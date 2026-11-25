#!/usr/bin/env python3
"""Per-patient methylation-change heterogeneity ("cladogram") metrics, independent re-implementation.

Per pair: probes -> units (fixed windows OR probe clusters), mean beta per unit in N and T,
delta = T - N.  Units are clustered in (beta_N, beta_T) space; the TRUNK cluster is the one with the
lowest variance of delta (clusters < --minclust units are ignored); the rest are PRIVATE.
Events = units with |delta| > 1 SD of delta within the sample (hyper if delta>0, hypo if <0).

  --mode windows : UPGMA on units (scipy), k by silhouette over 2..8 (or fixed --k)   [replica of thesis method]
  --mode clusters: probe clusters (gap <= --gap bp); Gaussian mixture (BIC, k 2..6) on units  [high resolution]

Outputs (--out): phylo_<mode>_patients.tsv, phylo_<mode>_chrom.tsv (events per chrom x trunk/private x hyper/hypo)
"""
import argparse, sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True); ap.add_argument("--betadir", required=True)
ap.add_argument("--manifest", required=True, help="HM450 manifest tsv(.gz): probe id, chrom, position columns")
ap.add_argument("--mode", choices=["windows", "clusters"], required=True)
ap.add_argument("--window", type=int, default=5_000_000); ap.add_argument("--gap", type=int, default=500)
ap.add_argument("--minprobes", type=int, default=None, help="default 20 (windows) / 2 (clusters)")
ap.add_argument("--minclust", type=int, default=10); ap.add_argument("--k", type=int, default=None)
ap.add_argument("--maxfit", type=int, default=100000, help="max units used to fit GMM")
ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--out", default="phylo_out")
ap.add_argument("--tcol", default="t_file"); ap.add_argument("--ncol", default="n_file")
ap.add_argument("--idcol", default="pid"); ap.add_argument("--typecol", default="type")
a = ap.parse_args()
minprobes = a.minprobes or (20 if a.mode == "windows" else 2)

# ---- manifest -> unit index per probe (done once)
hdr = pd.read_csv(a.manifest, sep="\t", nrows=0).columns
def pick(*c):
    low = {h.lower(): h for h in hdr}
    for x in c:
        if x.lower() in low: return low[x.lower()]
    sys.exit(f"manifest column not found among {c}; columns: {list(hdr)[:15]}")
cid, cch, cps = pick("probeID", "Probe_ID", "IlmnID"), pick("CpG_chrm", "chrom", "CHR"), pick("CpG_beg", "CpG_end", "MAPINFO", "pos")
cmask = next((h for h in hdr if h.lower() == "mask_general"), None)
cols = [cid, cch, cps] + ([cmask] if cmask else [])
m = pd.read_csv(a.manifest, sep="\t", usecols=cols, low_memory=False)[cols]
m.columns = ["id", "chr", "pos"] + (["mask"] if cmask else [])
m["chr"] = m.chr.astype(str).str.replace("^chr", "", regex=True)
m = m[m.chr.isin([str(i) for i in range(1, 23)] + ["X"]) & m.pos.notna()]
if cmask: m = m[~m["mask"].astype(str).str.upper().isin(["TRUE", "T", "1"])]
m["pos"] = m.pos.astype(int); m["ci"] = m.chr.map(lambda c: 23 if c == "X" else int(c))
m = m.sort_values(["ci", "pos"]).drop_duplicates("id").reset_index(drop=True)
if a.mode == "windows":
    key = m.ci * 10_000 + m.pos // a.window
    m["unit"] = pd.factorize(key)[0]
else:
    new = (m.ci.values[1:] != m.ci.values[:-1]) | (np.diff(m.pos.values) > a.gap)
    m["unit"] = np.concatenate([[0], np.cumsum(new)])
NU = int(m.unit.max()) + 1
unit_chr = m.groupby("unit").chr.first().values
ids, uidx = m.id.values, m.unit.values
print(f"{len(m)} probes after filtering -> {NU} units ({a.mode})", file=sys.stderr)

def read_beta(fn):
    p = Path(fn)
    p = p if p.is_absolute() else Path(a.betadir) / p
    d = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1], names=["id", "b"], dtype={"id": str}, low_memory=False)
    d = d[d.id.str.startswith("cg")]
    return pd.to_numeric(d.set_index("id")["b"], errors="coerce").reindex(ids).values

def unit_means(bn, bt):
    ok = ~(np.isnan(bn) | np.isnan(bt))
    cnt = np.bincount(uidx[ok], minlength=NU)
    sn = np.bincount(uidx[ok], weights=bn[ok], minlength=NU); st = np.bincount(uidx[ok], weights=bt[ok], minlength=NU)
    keep = cnt >= minprobes
    return sn[keep] / cnt[keep], st[keep] / cnt[keep], np.where(keep)[0], cnt[keep]

def cluster_windows(X):
    Z = linkage(X, "average")
    if a.k:
        return fcluster(Z, a.k, "maxclust"), a.k
    from sklearn.metrics import silhouette_score
    best = (-1, None, None)
    for k in range(2, 9):
        l = fcluster(Z, k, "maxclust")
        if len(set(l)) < 2: continue
        s = silhouette_score(X, l)
        if s > best[0]: best = (s, l, k)
    return best[1], best[2]

def cluster_gmm(X):
    from sklearn.mixture import GaussianMixture
    rng = np.random.default_rng(0)
    fit = X if len(X) <= a.maxfit else X[rng.choice(len(X), a.maxfit, replace=False)]
    ks = [a.k] if a.k else range(2, 7)
    best = None
    for k in ks:
        g = GaussianMixture(k, covariance_type="full", random_state=0, n_init=1).fit(fit)
        b = g.bic(fit)
        if best is None or b < best[0]: best = (b, g, k)
    return best[1].predict(X) + 1, best[2]

def one(row):
    try:
        bn, bt = read_beta(row[a.ncol]), read_beta(row[a.tcol])
        mn, mt, uid, cnt = unit_means(bn, bt)
        d = mt - mn
        if len(d) < 3 * a.minclust:
            print(f"SKIP {row[a.idcol]}: only {len(d)} units with >= {minprobes} probes", file=sys.stderr); return None, None
        X = np.column_stack([mn, mt])
        lab, k = (cluster_windows(X) if a.mode == "windows" else cluster_gmm(X))
        ok = [c for c in np.unique(lab) if (lab == c).sum() >= a.minclust]
        if not ok:
            print(f"SKIP {row[a.idcol]}: no cluster >= minclust", file=sys.stderr); return None, None
        trunk = min(ok, key=lambda c: np.var(d[lab == c]))
        sd = d.std()
        ev = np.abs(d) > sd
        is_tr = lab == trunk
        hyper = d > 0
        te, pe = (ev & is_tr).sum(), (ev & ~is_tr).sum()
        res = dict(pid=row[a.idcol], type=row[a.typecol], n_units=len(d), sd_delta=sd, mean_delta=d.mean(),
                   n_events=int(ev.sum()), frac_hyper=float((ev & hyper).sum() / max(ev.sum(), 1)), k=k,
                   trunk_unit_frac=float(is_tr.mean()), trunk_events=int(te), private_events=int(pe),
                   private_fraction=float(pe / max(te + pe, 1)),
                   private_hyper_frac=float(((ev & ~is_tr & hyper).sum()) / max(pe, 1)),
                   trunk_hyper_frac=float(((ev & is_tr & hyper).sum()) / max(te, 1)))
        ch = unit_chr[uid]
        rows = []
        for c in pd.unique(ch):
            s = ch == c
            for nm, mask in (("trunk", is_tr), ("private", ~is_tr)):
                rows.append(dict(pid=row[a.idcol], type=row[a.typecol], chrom=c, clone=nm,
                                 n_hyper=int((ev & mask & hyper & s).sum()), n_hypo=int((ev & mask & ~hyper & s).sum()),
                                 n_units=int((mask & s).sum())))
        return res, rows
    except Exception as e:
        print(f"FAILED {row[a.idcol]}: {e}", file=sys.stderr); return None, None

if __name__ == "__main__":
    pairs = pd.read_csv(a.pairs, sep="\t")
    recs = [r for _, r in pairs.iterrows()]
    with Pool(a.jobs) as pool:
        out = []
        for i, r in enumerate(pool.imap(one, recs, chunksize=1)):
            out.append(r)
            if (i + 1) % 25 == 0: print(f"{i+1}/{len(recs)}", file=sys.stderr)
    res = [o[0] for o in out if o[0]]; ch = [x for o in out for x in (o[1] or [])]
    Path(a.out).mkdir(parents=True, exist_ok=True)
    pd.DataFrame(res).to_csv(f"{a.out}/phylo_{a.mode}_patients.tsv", sep="\t", index=False)
    pd.DataFrame(ch).to_csv(f"{a.out}/phylo_{a.mode}_chrom.tsv", sep="\t", index=False)
    print(f"{len(res)}/{len(recs)} pairs ok -> {a.out}/")
