#!/usr/bin/env python3
"""C8.5 grade / stage / subtype availability per type for our 506 pairs.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c8_5_grade_subtype.py
Sources: TCGA-CDR (Liu 2018; grade, stage), GDC diagnoses (grade fallback), TCGASubtype.20170308 (mRNA subtype).
PRAD Gleason not used (C8.5 decision). Usable = >=2 levels with >=8 patients each.
Outputs -> results/c8_prolif/clinical_patients.tsv, clinical_availability.tsv
"""
import re, numpy as np, pandas as pd
D, O = "data/c8", "results/c8_prolif"
log = lambda *a: print(*a, flush=True)
TYPES13 = ["KIRC","BRCA","PRAD","HNSC","THCA","LIHC","LUSC","KIRP","UCEC","BLCA","LUAD","COAD","ESCA"]
MINN = 8

dx = pd.read_csv(f"{D}/gdc_diagnoses.tsv", sep="\t").set_index("patient")
P = pd.DataFrame(index=dx.index); P["project"] = dx.project.str.replace("TCGA-", "")
P = P[P.project.isin(TYPES13)]

def norm_grade(x):
    s = str(x).strip().upper()
    if s in ("HIGH GRADE", "HIGH"): return "High"
    if s in ("LOW GRADE", "LOW"): return "Low"
    m = re.search(r"G([1-4])", s)
    return f"G{m.group(1)}" if m else np.nan
def norm_stage(x):
    m = re.search(r"STAGE\s*(IV|III|II|I)", str(x).upper())
    return m.group(1) if m else np.nan

C = pd.read_excel(f"{D}/tcga_cdr_clinical.xlsx", sheet_name=0).set_index("bcr_patient_barcode")
P["grade"] = C.histological_grade.reindex(P.index).map(norm_grade)
fb = dx.tumor_grade.reindex(P.index).map(norm_grade)
P["grade_src"] = np.where(P.grade.notna(), "CDR", np.where(fb.notna(), "GDC", ""))
P["grade"] = P.grade.fillna(fb)
P["stage"] = C.ajcc_pathologic_tumor_stage.reindex(P.index).map(norm_stage)
P["stage_late"] = P.stage.map({"I": 0, "II": 0, "III": 1, "IV": 1})

S = pd.read_csv(f"{D}/tcga_subtypes.tsv", sep="\t")
S["patient"] = S.pan_samplesID.str[:12] if "pan_samplesID" in S else S.iloc[:, 0].str[:12]
S = S.drop_duplicates("patient").set_index("patient")
P["subtype_mRNA"] = S.Subtype_mRNA.reindex(P.index) if "Subtype_mRNA" in S else np.nan
P.to_csv(f"{O}/clinical_patients.tsv", sep="\t")

rows = []
for t, d in P.groupby("project"):
    for v in ("grade", "stage", "stage_late", "subtype_mRNA"):
        vc = d[v].dropna().value_counts()
        rows.append(dict(type=t, variable=v, n=int(vc.sum()),
                         levels="; ".join(f"{k}:{n}" for k, n in vc.sort_index().items()),
                         usable=bool((vc >= MINN).sum() >= 2)))
A = pd.DataFrame(rows)
A.to_csv(f"{O}/clinical_availability.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
log(A[A.variable != "stage_late"].to_string(index=False))
log("\nusable (>=2 levels with >=%d): " % MINN)
log(A[A.usable].groupby("variable").type.apply(lambda s: ", ".join(s)).to_string())
log(f"\ngrade source: {P.grade_src.value_counts().to_dict()}")
