#!/usr/bin/env python3
"""Supplementary tables for C4 (X-chromosome methylation events).
Run from Meth_TCGA_GSE213478/:  python3 scripts/c4_supp_tables.py
Output: results/c4_xchrom/Supp_Tables_C4.xlsx (one sheet per table: bold title row, note row, table).
Reads existing outputs only; no refitting. Tested strata need >=10 patients."""
import argparse, os
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

D = "results/c4_xchrom"
ap = argparse.ArgumentParser()
ap.add_argument("--bytype", default=f"{D}/C4_2_x_enrichment_by_type.csv")
ap.add_argument("--strata", default=f"{D}/C4_1_3_x_enrichment_strata.csv")
ap.add_argument("--patients", default=f"{D}/C4_1_3_x_patients_merged.csv")
ap.add_argument("--cna", default=f"{D}/x_cna_by_patient.tsv")
ap.add_argument("--units", default="5 Mb windows", help="unit label for titles")
ap.add_argument("--out", default=f"{D}/Supp_Tables_C4.xlsx")
ap.add_argument("--minpairs", type=int, default=10)
a = ap.parse_args()

def rd(p):
    if not os.path.exists(p):
        print(f"WARNING missing {p} -> sheet skipped"); return None
    return pd.read_csv(p, sep="\t" if p.endswith(".tsv") else ",")

bt, st, pt, cn = rd(a.bytype), rd(a.strata), rd(a.patients), rd(a.cna)
U = a.units
E = {"all": "all events", "hyper": "hypermethylation", "hypo": "hypomethylation"}
def tidy(df, lead):
    out = df[lead].copy()
    for k, lab in E.items():
        out[f"Enrichment, {lab}"] = df[f"{k}_enrich" if f"{k}_enrich" in df else f"{k}_enrich_pooled"]
        out[f"q, {lab}"] = df[f"{k}_q"]
    return out.round(4)
sheets = []

# S-C4.1 availability
if pt is not None:
    g = pt.groupby("type")
    t = pd.DataFrame({"Pairs (n)": g.size(), "Female (n)": g.sex.apply(lambda s: (s == "female").sum()),
                      "Male (n)": g.sex.apply(lambda s: (s == "male").sum()),
                      "ASCAT available (n)": g.x_rel_cn.apply(lambda s: s.notna().sum()),
                      "Female, X LOH ≥50% (n)": g.female_xloh.apply(lambda s: (s == "F_xLOH").sum()),
                      "Female, X copy-neutral (n)": g.apply(lambda d: (d.x_cn_neutral & (d.sex == "female")).sum())})
    t = t.sort_values("Pairs (n)", ascending=False).reset_index().rename(columns={"type": "Cancer type"})
    sheets.append(("S-C4.1", "Table S-C4.1. Sex and X copy-number data per cancer type",
                   "Sex from GDC demographic.sex_at_birth. Copy number from GDC ASCAT allele-specific segments (primary tumour). "
                   "Female X copy-neutral = X LOH <20% and X copy number 0.75–1.25× autosomal mean.", [t]))

# S-C4.2 by type
if bt is not None:
    b = bt.rename(columns={"type": "Cancer type", "n_pairs": "Pairs (n)", "all_frac_pat_X_gt_A": "Patients with X rate > autosomal rate (fraction)"})
    t = tidy(b, ["Cancer type", "Pairs (n)", "Patients with X rate > autosomal rate (fraction)"])
    t = t.sort_values("Pairs (n)", ascending=False)
    sheets.append(("S-C4.2", f"Table S-C4.2. X-chromosome enrichment of methylation events by cancer type ({U})",
                   "Enrichment = (X events / X tested units) ÷ (autosomal events / autosomal tested units), pooled over patients; 1 = no enrichment. "
                   f"Paired Wilcoxon (X vs autosomal rate per patient), BH across types; tested if ≥{a.minpairs} pairs.", [t]))

# S-C4.3 / S-C4.4 strata
if st is not None:
    s = st.rename(columns={"type": "Cancer type", "stratum": "Stratum", "n": "Patients (n)"})
    s["Patients (n)"] = s["Patients (n)"].astype(int)
    lead = ["Cancer type", "Stratum", "Patients (n)"]
    sex = tidy(s[s.family == "sex"], lead)
    sheets.append(("S-C4.3", f"Table S-C4.3. X-chromosome enrichment by sex ({U})",
                   f"As Table S-C4.2, stratified by sex. q only for strata with ≥{a.minpairs} patients (BH within stratum across types).", [sex]))
    lab = {"F_xLOH": "Female, X LOH ≥50%", "F_noXLOH": "Female, X LOH <50%", "female": "Female, X copy-neutral"}
    f1 = tidy(s[s.family == "female_xloh"], lead); f2 = tidy(s[(s.family == "cn_neutral_by_sex") & (s.Stratum == "female")], lead)
    f = pd.concat([f1, f2]); f["Stratum"] = f.Stratum.map(lab)
    sheets.append(("S-C4.4", f"Table S-C4.4. X-chromosome enrichment in females by X copy-number state ({U})",
                   "Tests whether enrichment is explained by X copy loss. Male copy-number strata omitted: ASCAT male X copy number "
                   "was not interpretable (Table S-C4.5).", [f.sort_values(["Cancer type", "Stratum"])]))

# S-C4.5 ASCAT X QC by sex
if cn is not None:
    q = cn.groupby("sex")[["x_rel_cn", "x_loh_frac", "auto_mean_cn"]].describe(percentiles=[.1, .25, .5, .75, .9])
    q = q.drop(columns=[c for c in q.columns if c[1] in ("mean", "std", "min", "max")]).round(3)
    q.columns = [f"{ {'x_rel_cn': 'X CN / autosomal CN', 'x_loh_frac': 'X LOH fraction', 'auto_mean_cn': 'Autosomal mean CN'}[v]}, {s}" for v, s in q.columns]
    sheets.append(("S-C4.5", "Table S-C4.5. ASCAT X-chromosome copy-number summaries by sex",
                   "Expected: females ~1× autosomal CN with low X LOH; males ~0.5× with X LOH ~1. Female values agree; male values "
                   "(median 0.23×, LOH 0.45) are not biologically plausible, so male copy-number strata are not reported.",
                   [q.reset_index().rename(columns={"sex": "Sex"})]))

wb = Workbook(); wb.remove(wb.active)
for name, title, note, dfs in sheets:
    ws = wb.create_sheet(name); ws["A1"] = title; ws["A1"].font = Font(bold=True); ws["A2"] = note
    r = 4
    for df in dfs:
        for j, c in enumerate(df.columns, 1):
            ws.cell(r, j, str(c)).font = Font(bold=True)
        for i, row in enumerate(df.itertuples(index=False), r + 1):
            for j, v in enumerate(row, 1):
                ws.cell(i, j, None if (isinstance(v, float) and np.isnan(v)) else v)
        r += len(df) + 3
    for j in range(1, ws.max_column + 1):
        ws.column_dimensions[get_column_letter(j)].width = 20
os.makedirs(os.path.dirname(a.out), exist_ok=True); wb.save(a.out)
print(f"wrote {a.out}: {', '.join(s[0] for s in sheets)}")
