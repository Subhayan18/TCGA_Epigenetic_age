#!/usr/bin/env python3
"""Builds the C1 tables read by c1_supp_tables.py. Run from Meth_TCGA_GSE213478/:
    python3 scripts/c1_tables.py [--out results/c1_clocks]
Inputs: results/c1_clocks/c1_out/mitotic_clocks_pairs.tsv, results/c1_clocks/pair_table_corrected.csv.
Outputs: Table_mitotic_clocks_tumour_vs_normal.csv, Table_normal_mitotic_clock_vs_age.csv, accel_stats_corrected.csv"""
import argparse, pandas as pd, numpy as np
from scipy.stats import wilcoxon, spearmanr
from statsmodels.stats.multitest import multipletests
ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results/c1_clocks")
ap.add_argument("--mitotic", default="results/c1_clocks/c1_out/mitotic_clocks_pairs.tsv")
ap.add_argument("--horvath", default="results/c1_clocks/pair_table_corrected.csv"); a = ap.parse_args()
m = pd.read_csv(a.mitotic, sep="\t"); h = pd.read_csv(a.horvath)
tc = "cancer_type" if "cancer_type" in m else "type"
keep = m[tc].value_counts(); keep = keep[keep >= 10].index; m = m[m[tc].isin(keep)]
CL = ["epiTOC1", "epiTOC2", "epiTOC2_irS", "HypoClock", "stemTOC"]
rows = []
for c in CL:
    for t, g in m.groupby(tc):
        d = g[c + "_T"] - g[c + "_N"]
        rows.append([c, t, len(g), g[c + "_N"].median(), g[c + "_T"].median(), d.median(), (d > 0).mean(), wilcoxon(d).pvalue])
r = pd.DataFrame(rows, columns="clock cancer_type n median_normal median_tumour median_delta_T_minus_N fraction_pairs_T_gt_N p_wilcoxon".split())
r["q_BH"] = r.groupby("clock").p_wilcoxon.transform(lambda p: multipletests(p, method="fdr_bh")[1])
r.to_csv(f"{a.out}/Table_mitotic_clocks_tumour_vs_normal.csv", index=False)
cl4 = ["epiTOC1", "epiTOC2", "HypoClock", "stemTOC"]; rows = []
for t, g in m.groupby(tc):
    row = {"cancer_type": t, "n": len(g), "age_SD_years": g.age.std(), "epiTOC2_rate_normal_median": g.epiTOC2_irS_N.median(),
           "epiTOC2_rate_tumour_median": g.epiTOC2_irS_T.median()}
    for c in cl4:
        rho, p = spearmanr(g.age, g[c + "_N"]); row[c + "_rho"] = rho; row[c + "_p"] = p
    rows.append(row)
v = pd.DataFrame(rows)
for c in cl4: v[c + "_q_BH"] = multipletests(v[c + "_p"], method="fdr_bh")[1]
v.drop(columns=[c + "_p" for c in cl4]).to_csv(f"{a.out}/Table_normal_mitotic_clock_vs_age.csv", index=False)
h = h[h.type.isin(keep)]
ac = h.groupby("type").apply(lambda g: pd.Series({"n": len(g), "median_normal_DNAmAge": g.DNAmAge_Normal.median(), "median_tumour_DNAmAge": g.DNAmAge_Tumour.median(),
     "median_accel_T_minus_N": g.accel.median(), "fraction_pairs_T_gt_N": (g.accel > 0).mean(), "normal_MAE_years": g.normal_err.abs().mean(),
     "p_wilcoxon": wilcoxon(g.accel).pvalue})).reset_index().rename(columns={"type": "cancer_type"})
ac["q_BH"] = multipletests(ac.p_wilcoxon, method="fdr_bh")[1]; ac.to_csv(f"{a.out}/accel_stats_corrected.csv", index=False)
print("wrote 3 tables to", a.out)
