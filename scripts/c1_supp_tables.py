#!/usr/bin/env python3
"""C1 supplementary tables -> Supp_Tables_C1.xlsx (bold title row + note row + captioned tables per sheet).
Reads existing result files only (no refitting). Run from Meth_TCGA_GSE213478/:
    python3 scripts/c1_supp_tables.py [--root .] [--out results/c1_clocks/Supp_Tables_C1.xlsx]
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

ap = argparse.ArgumentParser(); ap.add_argument("--root", default="."); ap.add_argument("--out", default="results/c1_clocks/Supp_Tables_C1.xlsx")
a = ap.parse_args(); B = Book(); R = a.root

defs = pd.DataFrame([
    ["epiTOC1", 385, "mean β of CpGs", "Yang 2016 Genome Biol 17:205"],
    ["epiTOC2", 163, "estimated cumulative stem-cell divisions (TNSC); rate = TNSC/age", "Teschendorff 2020 Genome Med 12:56"],
    ["HypoClock", 678, "1 − mean β of solo-WCGW CpGs", "Teschendorff 2020"],
    ["stemTOC", 371, "95th percentile of β", "Zhu 2024 Nat Commun"],
    ["Horvath", 353, "published linear model (intercept 0.6955), mean-imputed missing β", "Horvath 2013 Genome Biol 14:R115"]],
    columns=["Clock", "CpGs (n)", "Definition", "Reference"])
cov, s_cov = load(R, "*coverage*v2*", "c1_out_v2/*coverage*", "mitotic_clocks_coverage*")
B.sheet("S-C1.1", "Supplementary Table C1.1 - Clock definitions and per-sample CpG coverage",
        "Pairs: 531 rebuilt (01A + 11, HM450). Coverage = clock CpGs with a non-missing β per sample (v2 rerun); min/median/max over samples.",
        [("A. Clock definitions", defs), ("B. CpG coverage per clock (numeric columns)", summ(cov) if cov is not None else None)], ["c1_mitotic_clocks_v2.py", s_cov])

t, s = load(R, "Table_mitotic_clocks_tumour_vs_normal.csv")
B.sheet("S-C1.2", "Supplementary Table C1.2 - Tumour minus normal mitotic-clock difference by cancer type",
        "13 types with ≥10 matched pairs (506 pairs). Δ = tumour − normal; paired Wilcoxon signed-rank test; BH q across types within clock.",
        [("Per type and clock", t)], [s])

t, s = load(R, "Table_normal_mitotic_clock_vs_age.csv")
B.sheet("S-C1.3", "Supplementary Table C1.3 - Normal-tissue mitotic clocks versus chronological age",
        "Within-type Spearman correlation in matched normal samples; BH q across all tests. epiTOC2 rate = TNSC/age.",
        [("Per type and clock", t)], [s])

t, s = load(R, "accel_stats_corrected.csv")
im, s2 = load(R, "*imput*sens*.csv", "*imput*.csv", "*imput*.tsv")
if im is None:  # C0.8 results were not saved as a file; values from the C0.8 analysis
    im = pd.DataFrame([["Horvath", "KIRC", "+11.3 y (p 2e-17)", "+11.3 y (p 2e-17)"], ["Horvath", "LUAD", "+3.1 y (p 0.89)", "+2.5 y (p 0.95)"],
                       ["Mitotic (all 4 clocks)", "KIRC", "p 2e-24 (stemTOC example)", "p 3e-24"], ["Mitotic (all 4 clocks)", "LUAD", "p 6e-5", "p 1.2e-4 (minimum possible exact p at n=14)"]],
                      columns=["Clock", "Cancer type", "All pairs", "Without LUAD 44-6145 and KIRC B0-5402"]); s2 = "C0.8 (values transcribed)"
B.sheet("S-C1.4", "Supplementary Table C1.4 - Horvath age acceleration (corrected labels) and imputation sensitivity",
        "Acceleration = DNAmAge(tumour) − DNAmAge(normal); Wilcoxon signed-rank, BH q. Imputation: pairs with >28 missing Horvath CpGs dropped (LUAD 44-6145, KIRC B0-5402).",
        [("A. Horvath acceleration by type", t), ("B. Imputation sensitivity", im)], [s, s2])

h, s = load(R, "Table_C14_hetero_mixed_models_v1.csv")
if h is not None:
    mc = next((c for c in h.columns if c.lower() in ("metric", "outcome")), None)
    if mc: h = h[~h[mc].astype(str).str.contains("k_icl")]
B.sheet("S-C1.5", "Supplementary Table C1.5 - Per-patient Δ mitotic clock versus heterogeneity of methylation change",
        "Model: metric ~ Δclock (z) + signed mean Δβ + |mean Δβ| + log sd Δβ + CPE + age + sex + type; BH q over all tests. k_icl omitted (degenerate). 500 bp clusters (n=490) and 5 Mb windows (|Δβ|>0.05).",
        [("Coefficient of Δclock per metric and clock", h)], [s])

ped = [("A. Entity medians by array", "Table_pediatric_entity_medians_by_array.csv"),
       ("B. Clock versus age within entity and array", "Table_pediatric_clock_vs_age_within_entity.csv"),
       ("C. Relapse pairs", "Table_pediatric_relapse_pairs.csv"),
       ("D. Relapse tests (exact Wilcoxon)", "Table_pediatric_relapse_tests.csv")]
tabs, src = [], []
for cap, f in ped:
    d, s = load(R, f); tabs.append((cap, d)); src.append(s)
B.sheet("S-C1.6", "Supplementary Table C1.6 - Pediatric mitotic clocks",
        "Clocks on CpGs present in ≥95% of samples on 450K, EPIC and EPICv2; all comparisons within array (array is confounded with entity). Entities with ≥10 samples; BH q.",
        tabs, src)
B.save(a.out)
