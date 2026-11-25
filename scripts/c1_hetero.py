#!/usr/bin/env python3
"""C1.4 heterogeneity metrics of methylation change per patient (replaces trunk/private fraction).
Units (probe clusters --mode clusters --gap 500, or fixed windows --mode windows --window N); delta = mean beta T - N per unit.
Per pair, for thresholds t in 0.05/0.10/0.20 (absolute |delta|):
  frac_changed_t, discord_t (share of changed units opposite to the majority direction),
  rel = delta/(1-beta_N) if delta>0 else delta/beta_N (fraction of the possible change; removes baseline dependence),
  q95_rel_t (a fully clonal change gives rel ~ purity, so q95_rel/CPE ~ 1), subclonal_t (share of changed units with |rel| < 0.5*q95_rel; purity cancels),
  k_icl (1-D Gaussian mixture on rel of changed units at t=0.10, k=1..4, ICL, every component weight >= 5%).
Plus mean_delta, sd_delta, n_units. Output: hetero_<mode>_patients.tsv
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
ap.add_argument("--maxfit", type=int, default=20000, help="max units used to fit GMM")
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


TS = (0.05, 0.10, 0.20)
def k_icl(x):
    from sklearn.mixture import GaussianMixture
    rng = np.random.default_rng(0)
    x = x if len(x) <= a.maxfit else x[rng.choice(len(x), a.maxfit, replace=False)]
    x = x.reshape(-1, 1); best = (np.inf, 1)
    for k in range(1, 5):
        g = GaussianMixture(k, random_state=0).fit(x)
        if g.weights_.min() < 0.05: continue
        z = np.clip(g.predict_proba(x), 1e-12, 1)
        icl = g.bic(x) - 2 * np.sum(z * np.log(z))
        if icl < best[0]: best = (icl, k)
    return best[1]

def one(row):
    try:
        bn, bt = read_beta(row[a.ncol]), read_beta(row[a.tcol])
        mn, mt, uid, cnt = unit_means(bn, bt)
        d = mt - mn
        rel = np.where(d > 0, d / np.clip(1 - mn, 0.05, None), d / np.clip(mn, 0.05, None))
        if len(d) < 100:
            print(f"SKIP {row[a.idcol]}: only {len(d)} units", file=sys.stderr); return None
        r = dict(pid=row[a.idcol], type=row[a.typecol], n_units=len(d), mean_delta=d.mean(), sd_delta=d.std())
        for t in TS:
            sel = np.abs(d) > t; c = d[sel]; rc = rel[sel]; n = len(c); up = (c > 0).sum()
            r[f"frac_changed_{t}"] = n / len(d)
            r[f"discord_{t}"] = min(up, n - up) / n if n >= 50 else np.nan
            if n >= 50:
                ab = np.abs(rc); q95 = np.quantile(ab, .95)
                r[f"q95_rel_{t}"] = q95
                r[f"subclonal_{t}"] = float((ab < 0.5 * q95).mean())
        c = rel[np.abs(d) > 0.10]
        r["k_icl"] = k_icl(np.abs(c)) if len(c) >= 200 else np.nan
        for s_, cc in (("up", c[c > 0]), ("down", -c[c < 0])):
            r[f"k_icl_{s_}"] = k_icl(cc) if len(cc) >= 200 else np.nan
        return r
    except Exception as e:
        print(f"FAILED {row[a.idcol]}: {e}", file=sys.stderr); return None

if __name__ == "__main__":
    pairs = pd.read_csv(a.pairs, sep="\t")
    recs = [r for _, r in pairs.iterrows()]
    out = []
    with Pool(a.jobs) as pool:
        for i, r in enumerate(pool.imap(one, recs), 1):
            if r: out.append(r)
            if i % 25 == 0: print(f"{i}/{len(recs)}", file=sys.stderr, flush=True)
    o = Path(a.out); o.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(out).to_csv(o / f"hetero_{a.mode}_patients.tsv", sep="\t", index=False)
    print(f"{len(out)}/{len(recs)} pairs ok -> {o}/")
