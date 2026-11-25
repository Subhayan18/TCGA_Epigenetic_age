#!/usr/bin/env python3
"""C3.5: per-type summary of the 5 Mb replica (delta = T - N, verified labels),
to compare against the thesis Figure 4 / Results 4.3 numbers."""
import pandas as pd

c = pd.read_csv("phylo_out/phylo_windows_chrom.tsv", sep="\t", dtype={"chrom": str})
p = pd.read_csv("phylo_out/phylo_windows_patients.tsv", sep="\t")
c["chrom"] = c["chrom"].str.replace("chr", "", regex=False)
c["ev"] = c["n_hyper"] + c["n_hypo"]

g = c.groupby(["type", "clone"])[["n_hyper", "n_hypo", "ev"]].sum().unstack("clone")
out = pd.DataFrame({
    "n_pairs": p.groupby("type").size(),
    "median_k": p.groupby("type")["k"].median(),
    "private_share_pooled": g[("ev", "private")] / (g[("ev", "private")] + g[("ev", "trunk")]),
    "private_frac_median": p.groupby("type")["private_fraction"].median(),
    "hyper_trunk": g[("n_hyper", "trunk")] / g[("ev", "trunk")],
    "hyper_private": g[("n_hyper", "private")] / g[("ev", "private")],
})

# X share vs mean autosomal share (thesis definition: share of all events per type)
s = c.groupby(["type", "chrom"])["ev"].sum()
s = s / s.groupby("type").transform("sum")
auto = s[s.index.get_level_values("chrom").isin([str(i) for i in range(1, 23)])]
out["X_share"] = s.xs("X", level="chrom")
out["mean_auto_share"] = auto.groupby("type").mean()
out["X_gt_auto"] = out["X_share"] > out["mean_auto_share"]

pd.set_option("display.width", 200)
print(out.sort_values("n_pairs", ascending=False).round(3).to_string())
out.to_csv("phylo_out/C3_5_orientation_summary.csv")
