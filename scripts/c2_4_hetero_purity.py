#!/usr/bin/env python3
"""C2.4: are 500 bp heterogeneity metrics driven by purity/composition?
Per type (>=10 pairs): Spearman metric vs CPE, dEpi, dFib, dIC (BH per metric family);
pooled within-type: metric ~ CPE + dFib + dIC + type FE (HC3); CPE>=0.6 vs all median shift.
LIHC excluded from composition terms (EpiDISH fails there, C2.1)."""
import sys, re, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as fdc

het_f, frac_f, out = sys.argv[1], sys.argv[2], sys.argv[3]
pat = r"([A-Z0-9]{2}-[A-Z0-9]{4})"
h = pd.read_csv(het_f, sep=None, engine="python")
key = next(c for c in h.columns if h[c].astype(str).str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}").mean() > .95)
h["patient"] = h[key].astype(str).str.extract(pat)[0]
f = pd.read_csv(frac_f, sep="\t")
for c in ("Epi", "Fib", "IC"): f["d" + c] = f["T_" + c] - f["N_" + c]
f = f[["patient", "type", "purity_CPE", "dEpi", "dFib", "dIC"]]
drop = {"patient", "type", "cancer_type", "age", "sex", "purity_CPE", "CPE"}
mets = [c for c in h.select_dtypes("number").columns if c not in drop and h[c].nunique() > 2]
d = f.merge(h[["patient"] + mets], on="patient")
d = d[d.groupby("type")["patient"].transform("size") >= 10]
d.loc[d.type == "LIHC", ["dEpi", "dFib", "dIC"]] = np.nan
print(f"merged {len(d)} pairs on '{key}'; metrics: {mets}")

rows = []
for m in mets:
    for t, g in d.groupby("type"):
        r = {"metric": m, "type": t, "n": g[m].notna().sum()}
        for cov in ("purity_CPE", "dEpi", "dFib", "dIC"):
            x = g[[m, cov]].dropna()
            r[f"rho_{cov}"], r[f"p_{cov}"] = stats.spearmanr(x[m], x[cov]) if len(x) >= 8 and x[m].nunique() > 2 else (np.nan, np.nan)
        hi = g[g.purity_CPE >= .6][m].dropna()
        r["median_all"], r["median_CPE_ge0.6"], r["n_CPE_ge0.6"] = g[m].median(), hi.median(), len(hi)
        rows.append(r)
s = pd.DataFrame(rows)
for cov in ("purity_CPE", "dEpi", "dFib", "dIC"):
    p = s[f"p_{cov}"]; s[f"q_{cov}"] = np.nan
    s.loc[p.notna(), f"q_{cov}"] = fdc(p.dropna(), method="bh")
s.to_csv(f"{out}/C2_4_hetero_vs_purity_by_type.csv", index=False, float_format="%.4g")

pooled = []
for m in mets:
    x = d.dropna(subset=[m, "purity_CPE", "dFib", "dIC"])
    X = pd.get_dummies(x["type"], drop_first=True, dtype=float)
    Z = (x[["purity_CPE", "dFib", "dIC"]] - x[["purity_CPE", "dFib", "dIC"]].mean()) / x[["purity_CPE", "dFib", "dIC"]].std()
    X = np.column_stack([np.ones(len(x)), Z.values, X.values]); y = ((x[m] - x[m].mean()) / x[m].std()).values
    XtXi = np.linalg.pinv(X.T @ X); b = XtXi @ X.T @ y; e = y - X @ b
    hh = np.einsum("ij,jk,ik->i", X, XtXi, X); V = XtXi @ (X.T @ (X * ((e / (1 - hh)) ** 2)[:, None])) @ XtXi
    se = np.sqrt(np.diag(V)); p = 2 * stats.t.sf(np.abs(b / se), len(y) - X.shape[1])
    for i, cov in enumerate(("CPE", "dFib", "dIC"), 1):
        pooled.append({"metric": m, "term": cov, "beta_std": b[i], "se_HC3": se[i], "p": p[i], "n": len(y)})
pm = pd.DataFrame(pooled); pm["q_BH"] = fdc(pm["p"], method="bh")
pm.to_csv(f"{out}/C2_4_hetero_vs_purity_pooled.csv", index=False, float_format="%.4g")
print(pm.to_string(index=False))
sig = s[[c for c in s.columns if c.startswith("q_")]].lt(.05).sum()
print("per-type tests with q<0.05:", sig.to_dict())
