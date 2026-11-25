#!/usr/bin/env python3
"""C4.2: X-chromosome event enrichment normalised by tested units per chromosome.
Rate = events / tested units (5 Mb windows). Enrichment = X rate / autosomal rate.
Per patient: paired Wilcoxon of X rate vs autosomal rate within type, BH across types."""
import sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

f = sys.argv[1] if len(sys.argv) > 1 else "phylo_out/phylo_windows_chrom.tsv"
c = pd.read_csv(f, sep="\t", dtype={"chrom": str})
c["chrom"] = c["chrom"].str.replace("chr", "", regex=False)
c = c[c["chrom"].isin([str(i) for i in range(1, 23)] + ["X"])]
c["grp"] = np.where(c["chrom"] == "X", "X", "A")

# per patient x {X, A}: events (all / hyper / hypo) and units, summed over clones
pp = c.groupby(["pid", "type", "grp"])[["n_hyper", "n_hypo", "n_units"]].sum().unstack("grp")
out = pd.DataFrame(index=pp.index)
for k, cols in {"all": ["n_hyper", "n_hypo"], "hyper": ["n_hyper"], "hypo": ["n_hypo"]}.items():
    for g in ["X", "A"]:
        out[f"{k}_rate_{g}"] = sum(pp[(col, g)] for col in cols) / pp[("n_units", g)]
out["units_X"], out["units_A"] = pp[("n_units", "X")], pp[("n_units", "A")]
out = out.reset_index()
out.to_csv("phylo_out/C4_2_x_rates_patients.csv", index=False)  # merge sex here for C4.1

def bh(p):
    p = np.asarray(p, float); m = np.isfinite(p); q = np.full_like(p, np.nan)
    r = p[m]; o = np.argsort(r); s = r[o] * len(r) / np.arange(1, len(r) + 1)
    q[np.where(m)[0][o]] = np.minimum(1, np.minimum.accumulate(s[::-1])[::-1]); return q

rows = []
for t, d in out.groupby("type"):
    row = {"type": t, "n_pairs": len(d), "units_X_median": d["units_X"].median()}
    for k in ["all", "hyper", "hypo"]:
        ex, ea = d[f"{k}_rate_X"], d[f"{k}_rate_A"]
        # pooled enrichment: (sum X events / sum X units) / (sum A events / sum A units)
        row[f"{k}_enrich_pooled"] = (ex * d.units_X).sum() / d.units_X.sum() / ((ea * d.units_A).sum() / d.units_A.sum())
        row[f"{k}_frac_pat_X_gt_A"] = (ex > ea).mean()
        diff = ex - ea
        row[f"{k}_p"] = wilcoxon(diff).pvalue if len(d) >= 10 and (diff != 0).any() else np.nan
    rows.append(row)
res = pd.DataFrame(rows).set_index("type")
for k in ["all", "hyper", "hypo"]:
    res[f"{k}_q"] = bh(res[f"{k}_p"])
res = res.sort_values("n_pairs", ascending=False)

pd.set_option("display.width", 250)
cols = ["n_pairs", "units_X_median", "all_enrich_pooled", "all_frac_pat_X_gt_A", "all_q",
        "hyper_enrich_pooled", "hyper_q", "hypo_enrich_pooled", "hypo_q"]
print(res[cols].to_string(float_format=lambda v: f"{v:.3g}"))
res.to_csv("phylo_out/C4_2_x_enrichment_by_type.csv")
