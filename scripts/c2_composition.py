#!/usr/bin/env python3
"""C2: cell-type composition (EpiDISH RPC, Python port) + TCGA purity + composition-adjusted clocks.

RPC = EpiDISH::epidish(method="RPC"): MASS::rlm (Huber k=1.345, MAD scale, maxit=50, acc=1e-4)
of beta on reference centroids; negative coefs -> 0; rescale to sum 1 (Teschendorff 2017 BMC Bioinf 18:105).
Refs: centEpiFibIC.m (all types), centEpiFibFatIC.m (BRCA, primary for BRCA).
Purity: Aran 2015 Nat Commun 6:8971 (TCGAbiolinks Tumor.purity), joined on TCGA-XX-XXXX-01A.

Clock adjustment (paired, per type): delta_clock ~ 1 + delta_Fib + delta_IC [+ delta_Fat for BRCA]
(Epi omitted: fractions sum to 1). Intercept = T-N clock difference at equal composition. HC3 SEs, BH over types.
"""
import argparse, os, re, sys
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy import stats

K_HUBER = 1.345


def rlm_huber(X, y, maxit=50, acc=1e-4):
    """MASS::rlm(method='M', psi.huber, scale.est='MAD', test.vec='resid')."""
    w = np.ones_like(y)
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ coef
    for _ in range(maxit):
        old = resid
        scale = np.median(np.abs(resid)) / 0.6745
        if scale == 0:
            break
        u = np.abs(resid / scale)
        w = np.minimum(1.0, K_HUBER / np.maximum(u, 1e-300))
        sw = np.sqrt(w)
        coef = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        resid = y - X @ coef
        conv = np.sqrt(np.sum((old - resid) ** 2) / max(1e-20, np.sum(old ** 2)))
        if conv <= acc:
            break
    return coef


def rpc(beta, ref):
    """beta: Series indexed by CpG; ref: DataFrame CpG x celltype. Returns (fractions, n_cpg)."""
    common = ref.index.intersection(beta.index)
    b = beta.loc[common]
    ok = b.notna().values
    y = b.values[ok]
    R = ref.loc[common].values[ok]
    X = np.column_stack([np.ones(len(y)), R])
    c = rlm_huber(X, y)[1:]
    c[c < 0] = 0
    s = c.sum()
    c = c / s if s > 0 else np.full_like(c, np.nan)
    return pd.Series(c, index=ref.columns), int(ok.sum())


def read_beta(path):
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], na_values=["NA", "", "NaN"])
    if not np.issubdtype(d[1].dtype, np.number):  # header line present
        d = d.iloc[1:]
        d[1] = pd.to_numeric(d[1], errors="coerce")
    return pd.Series(d[1].values, index=d[0].values)


REFS = {}


def init(ref_dir):
    REFS["EpiFibIC"] = pd.read_csv(os.path.join(ref_dir, "centEpiFibIC.m.tsv"), sep="\t", index_col=0)
    REFS["EpiFibFatIC"] = pd.read_csv(os.path.join(ref_dir, "centEpiFibFatIC.m.tsv"), sep="\t", index_col=0)


def work(job):
    path, is_brca = job
    beta = read_beta(path)
    out = {}
    f, n = rpc(beta, REFS["EpiFibIC"])
    out.update({f"EFI_{k}": v for k, v in f.items()}); out["EFI_ncpg"] = n
    if is_brca:
        f, n = rpc(beta, REFS["EpiFibFatIC"])
        out.update({f"EFFI_{k}": v for k, v in f.items()}); out["EFFI_ncpg"] = n
    return out


def detect_cols(p):
    pat = re.compile(r"([A-Z0-9]{2}-[A-Z0-9]{4})")
    pt = ty = None
    for c in p.columns:
        v = p[c].astype(str)
        if pt is None and c not in ("t_file", "n_file") and v.str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}", regex=True).mean() > 0.95:
            pt = c
        if ty is None and v.str.fullmatch(r"(TCGA-)?[A-Z]{3,5}").mean() > 0.95 and v.nunique() < 40:
            ty = c
    if pt is None or ty is None:
        sys.exit(f"Cannot detect patient/type columns in pairs file: {list(p.columns)}")
    return pt, ty, pat


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); m = ~np.isnan(p)
    if m.sum():
        from scipy.stats import false_discovery_control
        q[m] = false_discovery_control(p[m], method="bh")
    return q


def ols_hc3(y, X):
    n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    S = X * ((e / (1 - h)) ** 2)[:, None]
    V = XtXi @ (X.T @ S) @ XtXi
    se = np.sqrt(np.diag(V))
    t = b / se
    p = 2 * stats.t.sf(np.abs(t), df=n - k)
    return b, se, p


def wilcox(x):
    x = x[~np.isnan(x)]
    return stats.wilcoxon(x).pvalue if len(x) >= 3 and np.any(x != 0) else np.nan


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--pairs", default="pairs_gdc.tsv")
    a.add_argument("--beta-dir", default="/home/chattopa/data_storage/TCGA_GDC_rebuilt")
    a.add_argument("--ref-dir", default="c2_ref")
    a.add_argument("--out", default="c2_out")
    a.add_argument("--min-pairs", type=int, default=10)
    a.add_argument("--threads", type=int, default=8)
    a.add_argument("--clock", action="append", default=[],
                   help="NAME=FILE:TUMOURCOL:NORMALCOL (omit cols to list columns and exit)")
    a = a.parse_args()
    os.makedirs(a.out, exist_ok=True)

    p = pd.read_csv(a.pairs, sep="\t")
    pt, ty, pat = detect_cols(p)
    p["patient"] = p[pt].astype(str).str.extract(pat)[0]
    p["type"] = p[ty].astype(str).str.replace("TCGA-", "", regex=False)
    print(f"pairs: {len(p)}; patient col={pt}; type col={ty}")

    # ---- C2.1 fractions ----
    fr_file = os.path.join(a.out, "fractions_pairs.tsv")
    if os.path.exists(fr_file):
        wide = pd.read_csv(fr_file, sep="\t")
        print("reusing", fr_file)
    else:
        jobs = []
        for _, r in p.iterrows():
            for side, col in (("T", "t_file"), ("N", "n_file")):
                jobs.append((os.path.join(a.beta_dir, str(r[col])), r["type"] == "BRCA"))
        with Pool(a.threads, initializer=init, initargs=(a.ref_dir,)) as pool:
            res = pool.map(work, jobs, chunksize=4)
        rows = []
        for i, (_, r) in enumerate(p.iterrows()):
            d = {"patient": r["patient"], "type": r["type"]}
            for j, side in enumerate(("T", "N")):
                d.update({f"{side}_{k}": v for k, v in res[2 * i + j].items()})
            rows.append(d)
        wide = pd.DataFrame(rows)
        # primary set: EpiFibFatIC for BRCA, EpiFibIC otherwise
        for s in ("T", "N"):
            for c in ("Epi", "Fib", "IC", "Fat"):
                efi = wide.get(f"{s}_EFI_{c}", pd.Series(np.nan, index=wide.index))
                effi = wide.get(f"{s}_EFFI_{c}", pd.Series(np.nan, index=wide.index))
                wide[f"{s}_{c}"] = np.where(wide["type"] == "BRCA", effi, efi if c != "Fat" else np.nan)
        # ---- C2.2 purity ----
        pur = pd.read_csv(os.path.join(a.ref_dir, "tcga_purity_aran2015.tsv"), sep="\t")
        pur = pur[pur["Sample.ID"].str.endswith("-01A")].copy()
        pur["patient"] = pur["Sample.ID"].str.extract(pat)[0]
        pur = pur.drop_duplicates("patient")[["patient", "CPE", "ABSOLUTE", "ESTIMATE", "LUMP", "IHC"]]
        wide = wide.merge(pur.add_prefix("purity_").rename(columns={"purity_patient": "patient"}), on="patient", how="left")
        wide.to_csv(fr_file, sep="\t", index=False)

    keep = wide.groupby("type")["patient"].transform("size") >= a.min_pairs
    w = wide[keep].copy()
    print(f"types >= {a.min_pairs} pairs: {w['type'].nunique()} / {len(w)} pairs")

    # ---- C2.5 supp table: fractions by type x tissue; paired T vs N ----
    rows = []
    for t, g in w.groupby("type"):
        for c in ("Epi", "Fib", "IC", "Fat"):
            if g[f"T_{c}"].isna().all():
                continue
            d = (g[f"T_{c}"] - g[f"N_{c}"]).values
            rows.append({"type": t, "n_pairs": len(g), "fraction": c,
                         "N_median": g[f"N_{c}"].median(), "N_IQR": f"{g[f'N_{c}'].quantile(.25):.3f}-{g[f'N_{c}'].quantile(.75):.3f}",
                         "T_median": g[f"T_{c}"].median(), "T_IQR": f"{g[f'T_{c}'].quantile(.25):.3f}-{g[f'T_{c}'].quantile(.75):.3f}",
                         "delta_median": np.median(d), "p_wilcoxon": wilcox(d)})
    s = pd.DataFrame(rows); s["q_BH"] = bh(s["p_wilcoxon"])
    s.to_csv(os.path.join(a.out, "Supp_composition_by_type.csv"), index=False, float_format="%.4g")

    # sanity: tumour Epi vs CPE purity (expect rho > 0)
    rows = []
    for t, g in w.groupby("type"):
        g2 = g.dropna(subset=["T_Epi", "purity_CPE"])
        if len(g2) >= 5:
            r = stats.spearmanr(g2["T_Epi"], g2["purity_CPE"])
            rows.append({"type": t, "n": len(g2), "rho_TEpi_CPE": r.statistic, "p": r.pvalue})
        else:
            rows.append({"type": t, "n": len(g2), "rho_TEpi_CPE": np.nan, "p": np.nan})
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "check_Epi_vs_CPE.csv"), index=False, float_format="%.4g")
    print(f"purity CPE available: {w['purity_CPE'].notna().sum()}/{len(w)} tumours; "
          f"missing by type: {w[w['purity_CPE'].isna()]['type'].value_counts().to_dict()}")

    # ---- C2.3 composition-adjusted clocks ----
    out = []
    for spec in a.clock:
        name, rest = spec.split("=", 1)
        parts = rest.split(":")
        f = pd.read_csv(parts[0], sep=None, engine="python")
        if len(parts) < 3:
            print(f"{name}: columns in {parts[0]}: {list(f.columns)}"); continue
        tcol, ncol = parts[1], parts[2]
        key = next(c for c in f.columns if f[c].astype(str).str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}", regex=True).mean() > 0.95)
        f["patient"] = f[key].astype(str).str.extract(pat)[0]
        m = w.merge(f[["patient", tcol, ncol]], on="patient", how="inner")
        m["dclock"] = m[tcol] - m[ncol]
        print(f"{name}: merged {len(m)}/{len(w)} pairs on '{key}'")
        for t, g in m.groupby("type"):
            covs = ["Fib", "IC"] + (["Fat"] if t == "BRCA" else [])
            D = np.column_stack([np.ones(len(g))] + [(g[f"T_{c}"] - g[f"N_{c}"]).values for c in covs])
            y = g["dclock"].values
            b, se, pv = ols_hc3(y, D)
            hi = g[g["purity_CPE"] >= 0.6]
            r = {"clock": name, "type": t, "n": len(g),
                 "unadj_median_delta": np.median(y), "unadj_p_wilcoxon": wilcox(y),
                 "adj_delta": b[0], "adj_se_HC3": se[0], "adj_p": pv[0],
                 "frac_pairs_dEpi_gt0.1": np.mean(np.abs(g["T_Epi"] - g["N_Epi"]) > 0.1),
                 "n_CPE_ge0.6": len(hi),
                 "CPE_ge0.6_median_delta": np.median(hi["dclock"]) if len(hi) else np.nan,
                 "CPE_ge0.6_p_wilcoxon": wilcox(hi["dclock"].values) if len(hi) else np.nan}
            for c, bb, pp in zip(covs, b[1:], pv[1:]):
                r[f"coef_d{c}"] = bb; r[f"p_d{c}"] = pp
            out.append(r)
    if out:
        o = pd.DataFrame(out)
        for c in ("unadj_p_wilcoxon", "adj_p", "CPE_ge0.6_p_wilcoxon"):
            o[c + "_BH"] = o.groupby("clock")[c].transform(lambda x: pd.Series(bh(x.values), index=x.index))
        o.to_csv(os.path.join(a.out, "Table_clock_composition_adjusted.csv"), index=False, float_format="%.4g")
        print(o[["clock", "type", "n", "unadj_median_delta", "adj_delta", "adj_p_BH", "CPE_ge0.6_median_delta"]].to_string(index=False))


if __name__ == "__main__":
    main()
