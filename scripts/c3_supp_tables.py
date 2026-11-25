#!/usr/bin/env python3
"""Supplementary tables for C3 (heterogeneity of methylation change, 500 bp vs 5 Mb).
Run from Meth_TCGA_GSE213478/:  python3 scripts/c3_supp_tables.py
Output: results/c3_phylo/Supp_Tables_C3.xlsx (one sheet per table: bold title row, note row, table).
Reads existing outputs only; no refitting. Main analysis = types with >=10 pairs."""
import argparse, os
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

ap = argparse.ArgumentParser()
ap.add_argument("--h500", default="results/c3_phylo/hetero_v2/hetero_clusters_patients.tsv")
ap.add_argument("--h500v1", default="results/c3_phylo/hetero/hetero_clusters_patients.tsv")
ap.add_argument("--h5m", default="results/c3_phylo/hetero_w5M/hetero_windows_patients.tsv")
ap.add_argument("--orient", default="results/c3_phylo/C3_5_orientation_summary.csv")
ap.add_argument("--out", default="results/c3_phylo/Supp_Tables_C3.xlsx")
ap.add_argument("--minpairs", type=int, default=10)
a = ap.parse_args()

def rd(p, **k):
    if not os.path.exists(p):
        print(f"WARNING missing {p} -> sheet skipped"); return None
    return pd.read_csv(p, sep="\t" if p.endswith(".tsv") else ",", **k)

def miqr(s):
    s = s.dropna()
    return "" if s.empty else f"{s.median():.3f} [{s.quantile(.25):.3f}–{s.quantile(.75):.3f}]"

h5, h1, hw, orient = rd(a.h500), rd(a.h500v1), rd(a.h5m), rd(a.orient)
keep = h5.type.value_counts().loc[lambda s: s >= a.minpairs].index
order = h5[h5.type.isin(keep)].groupby("type").size().sort_values(ascending=False).index
main = h5[h5.type.isin(keep)]
TS = ["0.05", "0.1", "0.2"]
sheets = []

# S-C3.1 units / coverage
t = main.groupby("type").agg(n=("pid", "size"), med=("n_units", "median"), mn=("n_units", "min")).reindex(order)
t.columns = ["Pairs (n)", "Units per pair, median (500 bp)", "Units per pair, min (500 bp)"]
if hw is not None:
    w = hw[hw.type.isin(keep)].groupby("type").n_units.median().reindex(order)
    t["Units per pair, median (5 Mb)"] = w
low = main.loc[main.n_units < 0.5 * main.n_units.median(), ["pid", "type", "n_units"]]
low.columns = ["Patient", "Type", "Units (500 bp)"]
sheets.append(("S-C3.1", "Table S-C3.1. Analysis units and coverage per cancer type",
               "Units = probe clusters (HM450 probes ≤500 bp apart, ≥2 probes; mean β per unit) or fixed 5 Mb windows (≥20 probes). "
               "Second block: pairs with <50% of the median unit count (kept; excluded in sensitivity analyses).",
               [t.reset_index().rename(columns={"type": "Cancer type"}), low]))

# S-C3.2 metrics by type, 500 bp
cols = {"sd_delta": "SD of Δβ", "mean_delta": "Mean Δβ"}
for x in TS:
    cols[f"frac_changed_{x}"] = f"Fraction changed, |Δβ|>{x}"
    cols[f"discord_{x}"] = f"Discordance, |Δβ|>{x}"
t = pd.DataFrame({v: main.groupby("type")[k].apply(miqr) for k, v in cols.items()}).reindex(order)
t.insert(0, "Pairs (n)", main.groupby("type").size().reindex(order))
sheets.append(("S-C3.2", "Table S-C3.2. Heterogeneity of methylation change by cancer type (500 bp probe clusters)",
               "Median [IQR] across pairs. Δβ = tumour − matched normal mean β per unit. Discordance = share of changed units "
               "opposite to the patient's majority direction (0 = all changes in one direction, 0.5 = balanced).",
               [t.reset_index().rename(columns={"type": "Cancer type"})]))

# S-C3.3 metric diagnostics (why discordance)
met = ["sd_delta", "mean_delta"] + [f"{m}_{x}" for m in ("frac_changed", "discord", "subclonal", "q95_rel") for x in TS]
met = [m for m in met if m in main]
def within_rho(df, x, y):
    r = df.groupby("type").apply(lambda g: g[x].corr(g[y], method="spearman") if g[[x, y]].dropna().shape[0] >= a.minpairs else np.nan)
    return r.median()
t = pd.DataFrame({"Metric": met,
                  "ρ vs SD of Δβ (pooled)": [main[m].corr(main.sd_delta, method="spearman") for m in met],
                  "ρ vs SD of Δβ (median within type)": [within_rho(main, m, "sd_delta") for m in met],
                  "ρ vs mean Δβ (pooled)": [main[m].corr(main.mean_delta, method="spearman") for m in met]}).round(3)
sheets.append(("S-C3.3", "Table S-C3.3. Dependence of candidate heterogeneity metrics on global methylation change (500 bp)",
               "Spearman ρ. Metrics strongly tracking SD of Δβ measure overall magnitude, not heterogeneity. "
               "Discordance is the only metric reported in the main text.", [t]))

# S-C3.4 500 bp vs 5 Mb
if hw is not None:
    m = main.merge(hw, on=["pid", "type"], suffixes=("_500", "_5M"))
    rows = []
    for ty in order:
        g = m[m.type == ty]; r = {"Cancer type": ty, "Pairs (n)": len(g)}
        for x in ["0.05", "0.1"]:
            a5, aw = f"discord_{x}_500", f"discord_{x}_5M"
            r[f"Usable at 5 Mb, |Δβ|>{x} (n)"] = int(g[aw].notna().sum())
            r[f"Discordance 500 bp, |Δβ|>{x}"] = miqr(g[a5]); r[f"Discordance 5 Mb, |Δβ|>{x}"] = miqr(g[aw])
            ok = g[[a5, aw]].dropna()
            r[f"ρ 500 bp vs 5 Mb, |Δβ|>{x}"] = round(ok[a5].corr(ok[aw], method="spearman"), 3) if len(ok) >= a.minpairs else np.nan
        rows.append(r)
    sheets.append(("S-C3.4", "Table S-C3.4. Resolution comparison: 500 bp probe clusters vs 5 Mb windows",
                   "Usable = ≥50 changed units in the pair. At 5 Mb few pairs reach the threshold; 500 bp is the primary resolution.",
                   [pd.DataFrame(rows)]))

# S-C3.5 mixture components
blocks = []
k2 = h5[["k_icl", "k_icl_up", "k_icl_down"]].apply(lambda s: s.value_counts(dropna=False)).fillna(0).astype(int)
k2.index = [("NA" if pd.isna(i) else int(i)) for i in k2.index]
k2.columns = ["k, |relative change|", "k, gains only", "k, losses only"]
blocks.append(k2.rename_axis("Components (k)").reset_index())
if h1 is not None:
    ct = pd.crosstab(h1.set_index("pid").k_icl.rename("v1"), h5.set_index("pid").k_icl.rename("v2"))
    ct.index.name = "k (signed, v1) \\ k (|rel|, v2)"
    blocks.append(ct.reset_index())
sheets.append(("S-C3.5", "Table S-C3.5. Number of mixture components in relative methylation change (500 bp)",
               "1-D Gaussian mixture on changed units (|Δβ|>0.10), k = 1–4 by ICL, each component ≥5%. Fitting on signed change (v1) "
               "splits gains from losses (k = 2); on |change| (v2) nearly all pairs have k = 1: discrete clonal tiers are not resolvable "
               "from bulk single-region tumour–normal pairs. Metric not used further.", blocks))

# S-C3.6 thesis replica orientation
if orient is not None:
    o = orient.copy()
    o["≥10 pairs"] = o.n_pairs >= a.minpairs
    o = o.rename(columns={"type": "Cancer type", "n_pairs": "Pairs (n)", "private_share_pooled": "Private share (pooled)",
                          "hyper_trunk": "Hypermethylated share, trunk", "hyper_private": "Hypermethylated share, private",
                          "X_share": "X share of events", "mean_auto_share": "Mean autosome share"})
    o = o.drop(columns=[c for c in ["median_k", "private_frac_median", "X_gt_auto"] if c in o]).round(3)
    sheets.append(("S-C3.6", "Table S-C3.6. Re-implementation of the original 5 Mb UPGMA analysis (direction check)",
                   "Δβ = tumour − normal on verified files. Reproduces the direction of hyper/hypomethylation reported previously; "
                   "the trunk/private split depends on the clustering choices and is not used.", [o]))

wb = Workbook(); wb.remove(wb.active)
for name, title, note, dfs in sheets:
    ws = wb.create_sheet(name); ws["A1"] = title; ws["A1"].font = Font(bold=True)
    ws["A2"] = note; ws["A2"].alignment = Alignment(wrap_text=False)
    r = 4
    for df in dfs:
        for j, c in enumerate(df.columns, 1):
            ws.cell(r, j, str(c)).font = Font(bold=True)
        for i, row in enumerate(df.itertuples(index=False), r + 1):
            for j, v in enumerate(row, 1):
                ws.cell(i, j, None if (isinstance(v, float) and np.isnan(v)) else (round(v, 4) if isinstance(v, float) else v))
        r += len(df) + 3
    for j in range(1, ws.max_column + 1):
        ws.column_dimensions[get_column_letter(j)].width = 22
os.makedirs(os.path.dirname(a.out), exist_ok=True); wb.save(a.out)
print(f"wrote {a.out}: {', '.join(s[0] for s in sheets)}")
