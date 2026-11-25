#!/usr/bin/env python3
"""Supplementary Table 23: clock ~ heterogeneity models (500 bp). Run from Meth_TCGA_GSE213478/:
    python3 scripts/c3_hetero_models.py
Output: results/c3_phylo/Table_C14_hetero_models.csv
Model per (metric, clock): metric ~ d_clock + mean_dbeta + rms_dbeta + sd_dbeta + CPE + age + male [+ C(type)], 13 types, HC3, BH over table.
std_beta: all continuous variables z-scored within type. beta_raw: same model on raw units. rms_dbeta = sqrt(sd^2 + mean^2)
(exact magnitude of dbeta from the stored summaries; replace by mean|dbeta| if that was the original term). ESCA has no CPE: dropped by listwise deletion."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
rd = lambda p: pd.read_csv(p, sep="," if p.endswith(".csv") else "\t", low_memory=False)
def pid(df):
    for c in df.columns:
        s = df[c].astype(str)
        if s.str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}").mean() > .9:
            df = df.copy(); df.index = "TCGA-" + s.str.extract(r"([A-Z0-9]{2}-[A-Z0-9]{4})")[0]; return df[~df.index.duplicated()]
M, H = pid(rd("c1_out/mitotic_clocks_pairs.tsv")), pid(rd("pair_table_corrected.csv"))
F, G = pid(rd("c2_out/fractions_pairs.tsv")), pid(rd("results/c3_phylo/hetero_v2/hetero_clusters_patients.tsv"))
CC = pid(rd("data/gdc_clinical/case_clinical.tsv"))
CLK = ["epiTOC1", "epiTOC2", "HypoClock", "stemTOC"]
D = pd.DataFrame(index=G.index)
for k in CLK: D["d_" + k] = M[k + "_delta"]
D["d_Horvath"] = H.DNAmAge_Tumour - H.DNAmAge_Normal
D["CPE"] = F.purity_CPE
D["age"] = pd.to_numeric(CC[[c for c in CC.columns if "age" in c.lower()][0]], errors="coerce")
D["male"] = (CC[[c for c in CC.columns if c.lower() in ("sex", "gender")][0]].astype(str).str.lower() == "male").astype(float)
D["mean_dbeta"], D["sd_dbeta"], D["type"] = G.mean_delta, G.sd_delta, G.type
D["rms_dbeta"] = np.sqrt(G.sd_delta ** 2 + G.mean_delta ** 2)
MET = [f"{m}_{t}" for m in ("discord", "subclonal") for t in ("0.05", "0.1", "0.2")]
MET = [m.replace(".", "p") for m in MET]
for m in MET: D[m] = G[m.replace("0p", "0.")]
TY = ["KIRC", "BRCA", "PRAD", "HNSC", "THCA", "LIHC", "LUSC", "KIRP", "UCEC", "BLCA", "LUAD", "COAD", "ESCA"]
D = D[D.type.isin(TY)]
cont = [c for c in D.columns if c not in ("type", "male")]
Z = D.copy()
for c in cont: Z[c] = D.groupby("type")[c].transform(lambda v: (v - v.mean()) / v.std(ddof=0) if v.std(ddof=0) > 0 else v * 0)
rows = []
for m in MET:
    for k in CLK + ["Horvath"]:
        f = f"{m} ~ d_{k} + mean_dbeta + rms_dbeta + sd_dbeta + CPE + age + male + C(type)"
        v = [m, f"d_{k}", "mean_dbeta", "rms_dbeta", "sd_dbeta", "CPE", "age", "male"]
        try:
            rz = smf.ols(f, Z.dropna(subset=v)).fit(cov_type="HC3"); rr = smf.ols(f, D.dropna(subset=v)).fit(cov_type="HC3")
            rows.append(dict(metric=m, clock=k, n=int(rz.nobs), std_beta=rz.params[f"d_{k}"], se=rz.bse[f"d_{k}"],
                             p=rz.pvalues[f"d_{k}"], beta_raw=rr.params[f"d_{k}"], se_raw=rr.bse[f"d_{k}"], r2=rz.rsquared))
        except Exception as e: print("failed", m, k, e)
T = pd.DataFrame(rows); T["q_BH"] = multipletests(T.p, method="fdr_bh")[1]
T.to_csv("results/c3_phylo/Table_C14_hetero_models.csv", index=False)
print(T.round(4).to_string(index=False))
