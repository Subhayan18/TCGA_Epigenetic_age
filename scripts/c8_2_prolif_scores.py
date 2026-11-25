#!/usr/bin/env python3
"""C8.2 proliferation scores per RNA sample. Run from Meth_TCGA_GSE213478/:
    python3 scripts/c8_2_prolif_scores.py
Needs pandas, numpy, scipy. Inputs from data/c8/ (c8_1_download.py). Outputs -> results/c8_prolif/
Scores (log2 TPM+1):
  core   mean z of MKI67, PCNA, MCM2-7 (z within cancer type, T and N together, so T-N deltas are on one scale)
  yang9  mean expression of the Yang 2016 mitotic index (CDKN3 ILF2 KDELR2 RFC4 TOP2A MCM3 KPNA2 CKS2 CDC2=CDK1)
  e2f / g2m  singscore-style: mean within-sample rank of set genes, scaled to [0,1] (Foroutan 2018, unidirectional)
"""
import json, os, numpy as np, pandas as pd
from scipy.stats import spearmanr

D, O = "data/c8", "results/c8_prolif"
os.makedirs(O, exist_ok=True)
log = lambda *a: print(*a, flush=True)

X = pd.read_csv(f"{D}/rna_log2tpm_proteincoding.tsv.gz", sep="\t", index_col=0)
sym = X.index.str.split("|").str[1]
X = X.groupby(sym).mean()                                   # collapse duplicate symbols
man = pd.read_csv(f"{D}/rna_manifest.tsv", sep="\t")
man["col"] = man.patient + "_" + man.side
meta = man.set_index("col")[["patient", "side", "project"]].loc[X.columns]
log(f"matrix {X.shape[0]} genes x {X.shape[1]} samples")

sets = {"core": ["MKI67", "PCNA"] + [f"MCM{i}" for i in range(2, 8)],
        "yang9": ["CDKN3", "ILF2", "KDELR2", "RFC4", "TOP2A", "MCM3", "KPNA2", "CKS2", "CDK1"]}
for k, f in [("e2f", "HALLMARK_E2F_TARGETS"), ("g2m", "HALLMARK_G2M_CHECKPOINT")]:
    sets[k] = list(json.load(open(f"{D}/{f}.json")).values())[0]["geneSymbols"]
for k, g in sets.items():
    miss = [x for x in g if x not in X.index]
    log(f"{k}: {len(g) - len(miss)}/{len(g)} genes present" + (f" (missing {miss[:10]})" if miss else ""))
    sets[k] = [x for x in g if x in X.index]

S = pd.DataFrame(index=X.columns)
Z = X.loc[sets["core"]].T.groupby(meta.project).transform(lambda v: (v - v.mean()) / v.std(ddof=0))
S["core"] = Z.mean(axis=1)
S["yang9"] = X.loc[sets["yang9"]].mean()
R = X.rank(axis=0)                                          # within-sample ranks
n = X.shape[0]
for k in ("e2f", "g2m"):
    m = len(sets[k])
    lo, hi = (m + 1) / 2, n - (m - 1) / 2                   # min/max possible mean rank
    S[k] = (R.loc[sets[k]].mean() - lo) / (hi - lo)
S = S.join(meta)
S.to_csv(f"{O}/prolif_scores_samples.tsv", sep="\t")

# wide per patient: tumour scores, normal scores, delta (complete pairs only)
T = S[S.side == "T"].set_index("patient"); N = S[S.side == "N"].set_index("patient")
sc = ["core", "yang9", "e2f", "g2m"]
W = T[["project"] + sc].add_suffix("_T").rename(columns={"project_T": "project"})
W = W.join(N[sc].add_suffix("_N"), how="left")
for s in sc: W[f"{s}_D"] = W[f"{s}_T"] - W[f"{s}_N"]
W.to_csv(f"{O}/prolif_scores_patients.tsv", sep="\t")
log(f"patients: {len(W)} with tumour RNA, {W.core_N.notna().sum()} with T+N")

# agreement among scores: within-type Spearman in tumours (median over types) and pooled
rows = []
for a in sc:
    for b in sc:
        if a < b:
            per = [spearmanr(g[f"{a}_T"], g[f"{b}_T"])[0] for _, g in W.groupby("project") if len(g) >= 10]
            rows.append(dict(a=a, b=b, median_within_type_rho=np.median(per), min_rho=np.min(per)))
A = pd.DataFrame(rows); A.to_csv(f"{O}/prolif_score_agreement.tsv", sep="\t", index=False)
log("agreement (tumours, within type):\n" + A.round(2).to_string(index=False))

# T > N check per type (expected for all; Yang 2016 Fig S6-S7)
chk = W.dropna(subset=["core_N"]).groupby("project").apply(
    lambda g: pd.Series({"n": len(g), **{f"{s}_T>N": (g[f"{s}_D"] > 0).mean() for s in sc}}))
chk.to_csv(f"{O}/prolif_T_gt_N.tsv", sep="\t")
log("fraction of pairs with tumour > normal:\n" + chk.round(2).to_string())
