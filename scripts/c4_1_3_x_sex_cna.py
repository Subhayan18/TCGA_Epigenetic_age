#!/usr/bin/env python3
"""C4.1 + C4.3: X enrichment (X event rate / autosomal rate, from C4.2) stratified by sex
and by X copy-number state. Strata tested if n >= 10; BH within each stratum across types.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c4_1_3_x_sex_cna.py [rates.csv]"""
import sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

D = "results/c4_xchrom"
r = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else f"{D}/C4_2_x_rates_patients.csv")
x = pd.read_csv(f"{D}/x_cna_by_patient.tsv", sep="\t")[["pid", "sex", "x_rel_cn", "x_loh_frac"]]
d = r.merge(x, on="pid", how="left")
print("sex missing:", d.sex.isna().sum(), "| ASCAT missing:", d.x_rel_cn.isna().sum())

# X copy-number neutral: female X ~ autosomal ploidy, no X LOH; male X ~ half ploidy
F, M = d.sex == "female", d.sex == "male"
d["x_cn_neutral"] = (F & d.x_loh_frac.lt(0.2) & d.x_rel_cn.between(0.75, 1.25)) | (M & d.x_rel_cn.between(0.375, 0.625))
d["female_xloh"] = np.where(F & d.x_loh_frac.ge(0.5), "F_xLOH", np.where(F & d.x_loh_frac.lt(0.5), "F_noXLOH", None))

def bh(p):
    p = np.asarray(p, float); m = np.isfinite(p); q = np.full_like(p, np.nan)
    v = p[m]; o = np.argsort(v); s = v[o] * len(v) / np.arange(1, len(v) + 1)
    q[np.where(m)[0][o]] = np.minimum(1, np.minimum.accumulate(s[::-1])[::-1]); return q

def summ(g):
    row = {"n": len(g)}
    for k in ["all", "hyper", "hypo"]:
        ex, ea = g[f"{k}_rate_X"], g[f"{k}_rate_A"]
        row[f"{k}_enrich"] = (ex * g.units_X).sum() / g.units_X.sum() / ((ea * g.units_A).sum() / g.units_A.sum())
        diff = ex - ea
        row[f"{k}_p"] = wilcoxon(diff).pvalue if len(g) >= 10 and (diff != 0).any() else np.nan
    return pd.Series(row)

strata = {"sex": d.sex, "female_xloh": d.female_xloh, "cn_neutral_by_sex": d.sex.where(d.x_cn_neutral)}
res = []
for name, lab in strata.items():
    t = d.assign(stratum=lab).dropna(subset=["stratum"]).groupby(["type", "stratum"]).apply(summ).reset_index()
    t.insert(0, "family", name)
    for k in ["all", "hyper", "hypo"]:
        t[f"{k}_q"] = np.nan
        for s in t.stratum.unique():
            i = t.stratum == s; t.loc[i, f"{k}_q"] = bh(t.loc[i, f"{k}_p"])
    res.append(t)
res = pd.concat(res)
res.to_csv(f"{D}/C4_1_3_x_enrichment_strata.csv", index=False)
d.to_csv(f"{D}/C4_1_3_x_patients_merged.csv", index=False)

pd.set_option("display.width", 250)
cols = ["family", "type", "stratum", "n", "all_enrich", "all_q", "hyper_enrich", "hyper_q", "hypo_enrich", "hypo_q"]
print(res[res.n >= 10][cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
