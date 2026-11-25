#!/usr/bin/env python3
"""C8.4 telomere length (Barthel 2017 ST1, 'Unpaired Set') for our pairs.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c8_4_telomere.py
Primary: log2(TL tumour / TL solid normal), same library type (+ centre recorded); tumour TL secondary.
Blood normals (NB) not used. Outputs -> results/c8_prolif/telomere_patients.tsv, telomere_coverage_by_type.tsv
"""
import re, numpy as np, pandas as pd
from scipy.stats import spearmanr
D, O = "data/c8", "results/c8_prolif"
log = lambda *a: print(*a, flush=True)
TYPES13 = ["KIRC","BRCA","PRAD","HNSC","THCA","LIHC","LUSC","KIRP","UCEC","BLCA","LUAD","COAD","ESCA"]
PRIO = {"WGS": 0, "LPS": 1, "WXS": 2}            # preferred library when several measurements exist

ids = {m[:12] for m in re.findall(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", open("data/pairs/pairs_gdc.tsv").read())}
proj = pd.read_csv(f"{D}/gdc_diagnoses.tsv", sep="\t").set_index("patient").project.str.replace("TCGA-", "")
U = pd.read_excel(f"{D}/barthel2017_telomere_ST1.xlsx", sheet_name="Unpaired Set")
U = U[U.PatientID.isin(ids) & U.SampleType.isin(["TP", "NT"])].copy()
U["code"] = U.SampleID.str[13:15]
U = U[((U.SampleType == "TP") & (U.code == "01")) | ((U.SampleType == "NT") & (U.code == "11"))]
U["prio"] = U.LibraryType.map(PRIO).fillna(9)
log(f"TL rows for our patients: {len(U)}; library types {U.LibraryType.value_counts().to_dict()}")

rows = []
for pid, g in U.groupby("PatientID"):
    T, N = g[g.SampleType == "TP"], g[g.SampleType == "NT"]
    r = dict(patient=pid)
    if len(T):
        t = T.sort_values("prio").iloc[0]
        r.update(tl_T=t.TL, logTL_T=np.log2(t.TL) if t.TL > 0 else np.nan, lib_T=t.LibraryType, center_T=t.Center, cov_T=t.Coverage)
    common = sorted(set(T.LibraryType) & set(N.LibraryType), key=lambda x: PRIO.get(x, 9))
    if common:                                          # ratio only within one library type
        lib = common[0]
        t = T[T.LibraryType == lib].iloc[0]; n = N[N.LibraryType == lib].iloc[0]
        r.update(tl_N=n.TL, lib_ratio=lib, center_ratio=t.Center, same_center=t.Center == n.Center,
                 log2_TN=np.log2(t.TL / n.TL) if t.TL > 0 and n.TL > 0 else np.nan)
    elif len(N):
        r.update(ratio_note="T and N on different library types")
    rows.append(r)
P = pd.DataFrame(rows).set_index("patient")
P["project"] = proj.reindex(P.index)
P = P[P.project.isin(TYPES13)]
W = pd.read_csv(f"{O}/prolif_scores_patients.tsv", sep="\t", index_col=0)
P = P.join(W[["yang9_T"]])
P.to_csv(f"{O}/telomere_patients.tsv", sep="\t")

g = P.groupby("project")
cov = pd.DataFrame({"tumour_TL": g.tl_T.count(), "ratio_T/N": g.log2_TN.count(),
                    "ratio_diff_centre": g.same_center.apply(lambda s: (s == False).sum()),
                    "median_log2_TN": g.log2_TN.median().round(2),
                    "frac_T_shorter": g.log2_TN.apply(lambda s: (s.dropna() < 0).mean()).round(2)}).reindex(TYPES13)
cov.to_csv(f"{O}/telomere_coverage_by_type.tsv", sep="\t")
log("coverage (expect tumour shorter than normal in most pairs; Barthel 2017):\n" + cov.to_string())
log(f"ratio library types: {P.lib_ratio.value_counts().to_dict()}; "
    f"T/N on different library types (no ratio): {(P.get('ratio_note').notna().sum() if 'ratio_note' in P else 0)}")
x = P[["log2_TN", "yang9_T"]].dropna()
log(f"pooled rho log2(T/N) vs tumour proliferation (yang9): {spearmanr(x.log2_TN, x.yang9_T)[0]:.2f} (n={len(x)})")
