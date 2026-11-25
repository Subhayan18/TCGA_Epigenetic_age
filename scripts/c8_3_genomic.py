#!/usr/bin/env python3
"""C8.3 genomic proxies per tumour: ABSOLUTE ploidy / WGD, aneuploidy score (Taylor 2018), PanCanAtlas QC flags.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c8_3_genomic.py
Matching: 15-char sample barcode of our methylation tumour (TCGA-XX-XXXX-01). Outputs -> results/c8_prolif/
"""
import re, numpy as np, pandas as pd
from scipy.stats import spearmanr
D, O = "data/c8", "results/c8_prolif"
log = lambda *a: print(*a, flush=True)
col = lambda df, pat: next(c for c in df.columns if re.search(pat, c, re.I))

# our tumours: patient -> methylation tumour sample (15 char); project from GDC diagnoses
bc = re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z]?")
tum = {}
for line in open("data/pairs/pairs_gdc.tsv"):
    h = bc.findall(line)
    if h: tum[h[0][:12]] = h[0][:15]
dx = pd.read_csv(f"{D}/gdc_diagnoses.tsv", sep="\t").set_index("patient")
P = pd.DataFrame({"sample15": pd.Series(tum)}); P.index.name = "patient"
P["project"] = dx.project.str.replace("TCGA-", "").reindex(P.index)
TYPES13 = ["KIRC","BRCA","PRAD","HNSC","THCA","LIHC","LUSC","KIRP","UCEC","BLCA","LUAD","COAD","ESCA"]
P = P[P.project.isin(TYPES13)]
log(f"{len(P)} tumours in 13 types")

# ABSOLUTE
A = pd.read_csv(f"{D}/absolute_purity_ploidy.tsv", sep="\t")
A["s15"] = A[col(A, r"^array$|^sample$")].str[:15]
A = A.drop_duplicates("s15").set_index("s15")
P["abs_call"] = A[col(A, "call status")].reindex(P.sample15).values
P["abs_purity"] = A[col(A, r"^purity$")].reindex(P.sample15).values
P["ploidy"] = A[col(A, r"^ploidy$")].reindex(P.sample15).values
P["doublings"] = A[col(A, "doubling")].reindex(P.sample15).values
P["wgd"] = (P.doublings >= 1).where(P.doublings.notna()).astype("float")
P.loc[~P.abs_call.astype(str).str.contains("call"), ["ploidy", "doublings", "wgd"]] = np.nan  # keep all ABSOLUTE call types

# P.loc[P.abs_call.ne("called") & P.abs_call.notna(), ["ploidy", "doublings", "wgd"]] = np.nan  # use called solutions only

# aneuploidy score
N = pd.read_csv(f"{D}/aneuploidy_arm_calls.tsv", sep="\t")
N["s15"] = N[col(N, r"^sample")].str[:15]
N = N.drop_duplicates("s15").set_index("s15")
P["aneuploidy"] = N[col(N, "aneuploidy")].reindex(P.sample15).values

# PanCanAtlas QC: any 'Do_not_use' flag on this tumour sample (any platform)
Q = pd.read_csv(f"{D}/sample_quality_annotations.tsv", sep="\t", low_memory=False)
Q["s15"] = Q[col(Q, "aliquot_barcode")].str[:15]
dnu = col(Q, "do_not_use")
bad = Q[Q[dnu].astype(str).str.lower().isin(["true", "1"])]
P["qc_do_not_use_platforms"] = bad.groupby("s15")[col(Q, "platform")].apply(lambda x: ",".join(sorted(set(x.dropna().astype(str)))) or "unspecified").reindex(P.sample15).values
# P["qc_do_not_use_platforms"] = bad.groupby("s15")[col(Q, "platform")].apply(lambda x: ",".join(sorted(set(x)))).reindex(P.sample15).values

# merge proliferation (C8.2) for a sanity check (Taylor 2018: aneuploidy ~ proliferation)
W = pd.read_csv(f"{O}/prolif_scores_patients.tsv", sep="\t", index_col=0)
P = P.join(W[["yang9_T", "g2m_T"]])
P.to_csv(f"{O}/genomic_patients.tsv", sep="\t")

# coverage / descriptives by type
g = P.groupby("project")
cov = pd.DataFrame({"n": g.size(), "ploidy_called": g.ploidy.count(), "aneuploidy": g.aneuploidy.count(),
                    "WGD_frac": g.wgd.mean().round(2), "ploidy_med": g.ploidy.median().round(2),
                    "aneup_med": g.aneuploidy.median(),
                    "qc_flagged": g.qc_do_not_use_platforms.count()}).reindex(TYPES13)
cov.to_csv(f"{O}/genomic_coverage_by_type.tsv", sep="\t")
log("coverage and descriptives:\n" + cov.to_string())
log(f"ABSOLUTE call status: {P.abs_call.value_counts(dropna=False).to_dict()}")
if P.qc_do_not_use_platforms.notna().any():
    log("QC-flagged tumours:\n" + P.loc[P.qc_do_not_use_platforms.notna(), ["project", "qc_do_not_use_platforms"]].to_string())

# sanity: within-type Spearman of aneuploidy / ploidy vs proliferation (yang9)
rows = []
for t, d in g:
    for v in ("aneuploidy", "ploidy"):
        x = d[[v, "yang9_T"]].dropna()
        if len(x) >= 10: rows.append(dict(type=t, var=v, n=len(x), rho=spearmanr(x[v], x.yang9_T)[0]))
S = pd.DataFrame(rows).pivot(index="type", columns="var", values="rho").round(2)
log("within-type rho with yang9 (expected mostly >0; Taylor 2018):\n" + S.to_string())
log(f"aneuploidy vs WGD (pooled point-biserial rho): {spearmanr(P.aneuploidy, P.wgd, nan_policy='omit')[0]:.2f}")
