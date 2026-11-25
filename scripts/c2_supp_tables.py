#!/usr/bin/env python3
"""C2 supplementary tables -> Supp_Tables_C2.xlsx (bold title row + note row + captioned tables per sheet).
Reads existing result files only (no refitting). Run from Meth_TCGA_GSE213478/:
    python3 scripts/c2_supp_tables.py [--root .] [--out results/c2_composition/Supp_Tables_C2.xlsx]
Needs pandas + openpyxl. Files are located by name under --root (archive/ and legacy_student/ are skipped)."""
import argparse
import re, sys
from pathlib import Path
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

SKIP = ("archive", "legacy_student", "NOT_FOR_MANUSCRIPT")
TOK = {"delta": "Δ", "q": "BH q", "p": "p", "rho": "Spearman ρ", "n": "n", "se": "SE", "beta": "β",
       "T": "tumour", "N": "normal", "pct": "%", "cpe": "CPE", "ic": "immune", "fib": "fibroblast",
       "epi": "epithelial", "BH": "", "gt": ">", "minus": "−", "wilcoxon": "Wilcoxon p", "fat": "adipocyte", "irs": "rate", "med": "median", "iqr": "IQR", "type": "Cancer type"}

EXACT = {"cancer_type": "Cancer type", "type": "Cancer type", "p_wilcoxon": "Wilcoxon p", "q_BH": "BH q"}

def nice(c):
    s = str(c)
    if s in EXACT: return EXACT[s]
    if " " in s: return s
    return " ".join(x for x in (TOK.get(t, TOK.get(t.lower(), t)) for t in s.split("_")) if x).strip()

def find(root, *pats):
    for pat in pats:
        hits = [p for p in Path(root).rglob(pat) if p.is_file() and p.stat().st_size > 0
                and not any(s in str(p) for s in SKIP)]
        if hits: return sorted(hits, key=lambda p: len(p.parts))[0]
    return None

def read(p):
    return pd.read_csv(p, sep=None, engine="python")

def summ(df, by=None):
    """min/median/max of numeric columns, optionally by group."""
    num = df.select_dtypes("number")
    if by is None:
        return num.agg(["min", "median", "max"]).T.reset_index().rename(columns={"index": "column"})
    g = pd.concat([df[by], num], axis=1).groupby(by).agg(["min", "median", "max"])
    g.columns = [f"{a} {b}" for a, b in g.columns]
    return g.reset_index()

class Book:
    def __init__(self): self.wb = Workbook(); self.wb.remove(self.wb.active); self.index = []
    def sheet(self, name, title, note, tables, sources):
        """tables = [(caption, DataFrame|None)]; None/empty -> flagged as missing."""
        ws = self.wb.create_sheet(name[:31])
        ws.append([title]); ws["A1"].font = Font(bold=True, size=12)
        ws.append([note]); ws["A2"].alignment = Alignment(wrap_text=False, vertical="top"); ws["A2"].font = Font(italic=True)
        for cap, df in tables:
            ws.append([]); ws.append([cap]); ws.cell(ws.max_row, 1).font = Font(bold=True)
            if df is None or len(df) == 0:
                ws.append(["SOURCE NOT FOUND OR EMPTY - rerun the producing step"]); continue
            df = df.copy(); df.columns = [nice(c) for c in df.columns]
            ws.append(list(df.columns)); h = ws.max_row
            for c in ws[h]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="E8E8E8"); c.alignment = Alignment(wrap_text=True, vertical="top")
            for row in df.itertuples(index=False):
                ws.append([None if (isinstance(v, float) and v != v) else v for v in row])
            for j, col in enumerate(df.columns, 1):
                small = re.search(r"\b(p|q)\b|BH|p value|fdr", col, re.I)
                for r in range(h + 1, ws.max_row + 1):
                    c = ws.cell(r, j)
                    if isinstance(c.value, float):
                        c.number_format = "0.0E+00" if small else ("0.00" if abs(c.value) >= 10 else "0.000")
        for j in range(1, ws.max_column + 1):
            ws.column_dimensions[get_column_letter(j)].width = 16
        ws.column_dimensions["A"].width = 22; ws.freeze_panes = "A3"
        self.index.append((name, title, "; ".join(f"{s}" for s in sources)))
        print(f"{name}: " + ", ".join(f"{c.split(':')[0]}={'OK' if d is not None and len(d) else 'MISSING'}" for c, d in tables))
    def save(self, out):
        ws = self.wb.create_sheet("Index", 0); ws.append(["Sheet", "Content", "Source files"])
        for c in ws[1]: c.font = Font(bold=True)
        for r in self.index: ws.append(list(r))
        for col, w in zip("ABC", (12, 80, 90)): ws.column_dimensions[col].width = w
        Path(out).parent.mkdir(parents=True, exist_ok=True); self.wb.save(out); print("wrote", out)

def load(root, *pats):
    p = find(root, *pats)
    if p is None: print("  not found:", pats[0]); return None, pats[0]
    return read(p), p.name

ap = argparse.ArgumentParser(); ap.add_argument("--root", default="."); ap.add_argument("--out", default="results/c2_composition/Supp_Tables_C2.xlsx")
a = ap.parse_args(); B = Book(); R = a.root

fr, s_fr = load(R, "fractions_pairs.tsv")
tc = next((c for c in fr.columns if c.lower() in ("type", "cancer", "project")), None) if fr is not None else None
qc = av = None
if fr is not None and tc:
    cc = [c for c in fr.columns if re.search(r"cpg|nref|n_ref", c, re.I) and fr[c].dtype.kind in "if"]
    if cc:
        qc = summ(fr[[tc] + cc], by=tc)
    pc = next((c for c in fr.columns if re.search(r"cpe|purity", c, re.I)), None)
    if pc:
        g = fr.groupby(tc)[pc]
        av = pd.DataFrame({"pairs (n)": g.size(), "with CPE (n)": g.count(), "median CPE": g.median(),
                           "CPE ≥0.6 (n)": g.apply(lambda x: (x >= 0.6).sum())}).reset_index()
chk, s_chk = load(R, "check_Epi_vs_CPE.csv")
B.sheet("S-C2.1", "Supplementary Table C2.1 - Deconvolution quality control and purity availability",
        "EpiDISH RPC (centEpiFibIC.m, 716 CpGs; BRCA centEpiFibFatIC.m, 491 CpGs; Teschendorff 2017 BMC Bioinformatics 18:105). CPE = consensus purity estimate (Aran 2015 Nat Commun 6:8971); absent for ESCA and 4 BRCA. LIHC deconvolution fails validation (tumour epithelial fraction vs CPE ρ 0.13).",
        [("A. Reference CpGs used per sample (by type)", qc), ("B. CPE availability per type", av),
         ("C. Tumour epithelial fraction versus CPE (validation)", chk)], [s_fr, s_chk])

t, s = load(R, "Supp_composition_by_type.csv")
B.sheet("S-C2.2", "Supplementary Table C2.2 - Cell-type composition by cancer type and tissue",
        "Median and IQR of EpiDISH fractions in tumour and matched normal; paired Wilcoxon signed-rank, BH q. LIHC fractions not interpretable (hepatocyte reference absent).",
        [("Composition by type", t)], [s])

t, s = load(R, "Table_clock_composition_adjusted.csv")
B.sheet("S-C2.3", "Supplementary Table C2.3 - Composition-adjusted clock differences",
        "Δclock ~ 1 + ΔFibroblast + ΔImmune (+ ΔAdipocyte for BRCA); HC3 SE; BH q. Intercept = tumour − normal difference at zero compositional change; |Δepithelial|>0.1 in 74-94% of pairs in several types, so the intercept is an extrapolation. Sensitivity: CPE ≥0.6. LIHC adjusted estimates invalid.",
        [("Adjusted and CPE ≥0.6 estimates per type and clock", t)], [s])

h, s = load(R, "C2_4_hetero_vs_purity_pooled.csv")
bt, s2 = load(R, "C2_4_hetero_vs_purity_by_type.csv")
for d in (h, bt):
    if d is not None:
        mc = next((c for c in d.columns if c.lower() in ("metric", "outcome")), None)
        if mc: d.drop(d[d[mc].astype(str).str.contains("k_icl")].index, inplace=True)
B.sheet("S-C2.4", "Supplementary Table C2.4 - Heterogeneity metrics versus tumour purity and composition",
        "500 bp units, 459 pairs (12 types; ESCA no CPE; LIHC composition terms excluded). Pooled: metric ~ CPE + ΔFibroblast + ΔImmune + type; standardised; HC3; BH q. Discordance is independent of purity and composition; magnitude metrics scale with purity; subclonal share tracks immune infiltration (not interpreted). k_icl omitted (degenerate).",
        [("A. Pooled model", h), ("B. Within-type Spearman", bt)], [s, s2])
B.save(a.out)
