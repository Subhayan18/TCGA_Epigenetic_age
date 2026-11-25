#!/usr/bin/env python3
"""C8.10 outputs: Figure 6 (composite, NEJM style) + Supp_Tables_C8.xlsx. No model refitting.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c8_10_outputs.py
Inputs: results/c8_prolif/*.tsv and results/c8_prolif/models/*.tsv; c1_out/mitotic_clocks_pairs.tsv (type order)
Outputs: results/c8_prolif/figures/Figure6_C8.svg/.pdf, results/c8_prolif/Supp_Tables_C8.xlsx
Style (C5 plan): Helvetica/Arial 6.5-7 pt, editable SVG text, 7 in wide, NEJM palette,
fixed 13-type order (ascending median stemTOC delta), filled marker = BH q<0.05, open = ns.
"""
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

R = "results/c8_prolif/"; M = R + "models/"; FO = R + "figures/"; os.makedirs(FO, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5})
NEJM = {"red": "#BC3C29", "blue": "#0072B5", "orange": "#E18727", "green": "#20854E", "purple": "#7876B1", "grey": "#6F99AD"}
CLOCKS = ["epiTOC1", "epiTOC2", "HypoClock", "stemTOC", "Horvath"]
CCOL = {"epiTOC1": NEJM["red"], "epiTOC2": NEJM["orange"], "HypoClock": NEJM["blue"], "stemTOC": NEJM["green"], "Horvath": NEJM["purple"]}
rd = lambda f: pd.read_csv(f, sep="\t")

# fixed type order: ascending median stemTOC delta
mc = pd.read_csv("c1_out/mitotic_clocks_pairs.tsv", sep="\t")
idc = next(c for c in mc.columns if mc[c].astype(str).str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}").mean() > 0.9)
mc["patient"] = "TCGA-" + mc[idc].astype(str).str.extract(r"([A-Z0-9]{2}-[A-Z0-9]{4})")[0]
pr = rd(R + "prolif_scores_patients.tsv")[["patient", "project"]]
mc = mc.merge(pr, on="patient"); mc["d"] = mc.stemTOC_T - mc.stemTOC_N
TYPES = mc.groupby("project").d.median().sort_values().index.tolist()

P6, W6, Y6 = rd(M + "C8_6_pooled.tsv"), rd(M + "C8_6_within_type.tsv"), rd(M + "C8_6_supp_tumour_only_yang.tsv")
B8, I7 = rd(M + "C8_8_block_partialR2.tsv"), rd(M + "C8_7_incremental.tsv")

# ---------------- Figure 6 ----------------
fig = plt.figure(figsize=(7, 6.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1.05], wspace=0.62, hspace=0.75,
                      left=0.10, right=0.98, top=0.95, bottom=0.13)
lab = lambda ax, s: ax.text(-0.18, 1.04, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")

# A: forest, proliferation effect on dclock
ax = fig.add_subplot(gs[0, 0]); lab(ax, "A")
MODELS = [("pooled_main_yang9", "Main", "o"), ("sens_excl_KIRC", "Excluding KIRC", "s"),
          ("pooled_dRNA_yang9", "Δ proliferation (T+N RNA)", "^"), ("sens_CPE>=0.6", "CPE ≥ 0.6", "D")]
for i, k in enumerate(CLOCKS):
    for j, (m, name, mk) in enumerate(MODELS):
        r = P6[(P6.model == m) & (P6.clock == k)]
        if r.empty: continue
        r = r.iloc[0]; y = i + (j - 1.5) * 0.18
        ax.plot([r.beta - 1.96 * r.se, r.beta + 1.96 * r.se], [y, y], color=CCOL[k], lw=0.7)
        ax.plot(r.beta, y, mk, ms=3.2, mec=CCOL[k], mew=0.7, mfc=CCOL[k] if r.q < 0.05 else "white")
ax.axvline(0, color="0.6", lw=0.5, ls="--")
ax.set_yticks(range(len(CLOCKS))); ax.set_yticklabels(CLOCKS); ax.invert_yaxis()
ax.set_xlabel("β (SD Δclock per SD tumour proliferation), 95% CI")
ax.set_title("Proliferation effect on tumour–normal clock change", loc="left")
for j, (m, name, mk) in enumerate(MODELS):
    ax.plot([], [], mk, ms=3.2, mec="0.3", mfc="0.3", label=name)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, handletextpad=0.3, columnspacing=1)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# B: partial R2 heatmap (no-telomere scope, n=465)
ax = fig.add_subplot(gs[0, 1]); lab(ax, "B")
b = B8[B8.scope == "no_telomere"]
BL = ["proliferation", "composition", "global_dbeta", "genomic"]
BLN = {"proliferation": "Proliferation", "composition": "Composition", "global_dbeta": "Global Δβ", "genomic": "Genomic"}
H = b.pivot(index="block", columns="clock", values="partial_R2").reindex(index=BL, columns=CLOCKS) * 100
Q = b.pivot(index="block", columns="clock", values="q").reindex(index=BL, columns=CLOCKS)
cm = LinearSegmentedColormap.from_list("nejm_blue", ["#FFFFFF", NEJM["blue"]])
ax.imshow(H.values, cmap=cm, vmin=0, vmax=max(30, np.nanmax(H.values)), aspect="auto")
for i in range(len(BL)):
    for j in range(len(CLOCKS)):
        v = H.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.1f}%" + ("*" if Q.values[i, j] < 0.05 else ""), ha="center", va="center",
                    fontsize=6, color="white" if v > 18 else "black")
ax.set_xticks(range(len(CLOCKS))); ax.set_xticklabels(CLOCKS, rotation=30, ha="right")
ax.set_yticks(range(len(BL))); ax.set_yticklabels([BLN[x] for x in BL])
ax.set_title(f"Partial R² of Δclock by block (n={int(b.n.max())}; * q<0.05)", loc="left")
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)

# C: within-type proliferation effect
ax = fig.add_subplot(gs[1, 0]); lab(ax, "C")
Wb = W6.pivot(index="type", columns="clock", values="beta").reindex(index=TYPES, columns=CLOCKS)
Wq = W6.pivot(index="type", columns="clock", values="q").reindex(index=TYPES, columns=CLOCKS)
cmd = LinearSegmentedColormap.from_list("nejm_div", [NEJM["blue"], "#FFFFFF", NEJM["red"]])
lim = np.nanmax(np.abs(Wb.values)); lim = 0.7 if not np.isfinite(lim) else max(0.5, min(lim, 1.0))
im = ax.imshow(Wb.T.values, cmap=cmd, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for i in range(len(CLOCKS)):
    for j in range(len(TYPES)):
        if Wq.T.values[i, j] < 0.05: ax.text(j, i, "*", ha="center", va="center", fontsize=8)
ax.set_xticks(range(len(TYPES))); ax.set_xticklabels(TYPES, rotation=45, ha="right")
ax.set_yticks(range(len(CLOCKS))); ax.set_yticklabels(CLOCKS)
ax.set_title("Proliferation effect within cancer type (β; * q<0.05)", loc="left"); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.ax.tick_params(labelsize=5.5); cb.outline.set_linewidth(0.4)

# D: incremental value beyond proliferation + global dbeta
ax = fig.add_subplot(gs[1, 1]); lab(ax, "D")
OUT = [("KIRC grade (ordinal)", "KIRC grade", False), ("late stage pooled", "Late stage (pooled)", False),
       ("late stage pooled excl KIRC", "Late stage, excl. KIRC", False), ("WGD pooled", "WGD (pooled)", False),
       ("aneuploidy pooled (excl THCA, PRAD)", "Aneuploidy (pooled)", False),
       ("HNSC G3 vs G1-2", "HNSC grade†", True), ("LIHC G3 vs G1-2", "LIHC grade†", True)]
OUT = [o for o in OUT if o[0] in set(I7.outcome)]
XL = 2.5
for i, (o, name, grey) in enumerate(OUT):
    for j, k in enumerate(CLOCKS):
        r = I7[(I7.outcome == o) & (I7.clock == k)]
        if r.empty or not np.isfinite(r.beta_dclock.iloc[0]): continue
        r = r.iloc[0]; y = i + (j - 2) * 0.13; x = float(np.clip(r.beta_dclock, -XL, XL))
        c = "0.75" if grey else CCOL[k]
        ax.plot(x, y, ">" if r.beta_dclock > XL else ("<" if r.beta_dclock < -XL else "o"), ms=3.2, mec=c, mew=0.7,
                mfc=c if (r.q < 0.05 and not grey) else "white")
ax.axvline(0, color="0.6", lw=0.5, ls="--"); ax.set_xlim(-XL - 0.1, XL + 0.1)
ax.set_yticks(range(len(OUT))); ax.set_yticklabels([o[1] for o in OUT]); ax.invert_yaxis()
for t, o in zip(ax.get_yticklabels(), OUT):
    if o[2]: t.set_color("0.6")
ax.set_xlabel("Δclock effect beyond proliferation and global Δβ\n(per SD; arrows = clipped; † dropped)")
ax.set_title("Added information of clock change", loc="left")
for k in CLOCKS: ax.plot([], [], "o", ms=3.2, color=CCOL[k], label=k)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=3, handletextpad=0.1, columnspacing=0.6)
for s in ("top", "right"): ax.spines[s].set_visible(False)

fig.savefig(FO + "Figure6_C8.svg"); fig.savefig(FO + "Figure6_C8.pdf")
print(f"figure -> {FO}Figure6_C8.svg/.pdf; type order: {', '.join(TYPES)}")

# ---------------- Supplementary tables ----------------
def tab(sheet, title, note, df):
    return (sheet, title, note, df)
cov = []
pp = rd(R + "prolif_scores_patients.tsv")
ge, tl, cl = rd(R + "genomic_patients.tsv"), rd(R + "telomere_patients.tsv"), rd(R + "clinical_patients.tsv")
for t in TYPES:
    p = pp[pp.project == t]; g = ge[ge.project == t]; L = tl[tl.project == t]; c = cl[cl.project == t]
    cov.append({"Cancer type": t, "Tumour RNA (n)": len(p), "Tumour + normal RNA (n)": int(p.yang9_N.notna().sum()) if "yang9_N" in p else int(p.yang9_D.notna().sum()),
                "Ploidy/WGD (n)": int(g.ploidy.notna().sum()), "Aneuploidy score (n)": int(g.aneuploidy.notna().sum()),
                "Telomere T/N ratio (n)": int(L.log2_TN.notna().sum()), "Grade (n)": int(c.grade.notna().sum()),
                "Stage (n)": int(c.stage.notna().sum()) if "stage" in c else np.nan, "mRNA subtype (n)": int(c.subtype_mRNA.notna().sum())})
S1 = pd.DataFrame(cov)

opt = lambda f: rd(f) if os.path.exists(f) else pd.DataFrame({"note": [f"{f} not found"]})
ag = opt(R + "prolif_score_agreement.tsv").rename(columns={"a": "Score A", "b": "Score B", "median_within_type_rho": "Median within-type ρ (tumours)", "min_rho": "Minimum within-type ρ"})
tn = opt(R + "prolif_T_gt_N.tsv").rename(columns={"project": "Cancer type", "n": "Pairs with T+N RNA (n)"})
tn.columns = [c.replace("_T>N", ": fraction tumour > normal") for c in tn.columns]
defs = pd.DataFrame({"Score": ["yang9 (primary)", "g2m (sensitivity)", "e2f", "core"],
                     "Definition": ["Mean log2(TPM+1) of CDKN3, ILF2, KDELR2, RFC4, TOP2A, MCM3, KPNA2, CKS2, CDK1 (Yang 2016 mitotic index)",
                                    "Rank-based single-sample score, Hallmark G2M_CHECKPOINT (singscore-style)",
                                    "Rank-based single-sample score, Hallmark E2F_TARGETS (not reported separately; ρ 0.97 with g2m)",
                                    "Mean within-type z of MKI67, PCNA, MCM2–7 (supplementary only)"]})

gc = rd(R + "genomic_coverage_by_type.tsv").rename(columns={"project": "Cancer type", "n": "Pairs (n)", "ploidy_called": "Ploidy called (n)",
     "aneuploidy": "Aneuploidy score (n)", "WGD_frac": "WGD fraction", "ploidy_med": "Median ploidy", "aneup_med": "Median aneuploidy score", "qc_flagged": "PanCanAtlas QC-flagged (n)"})
qcf = pd.DataFrame({"Patient": ["TCGA-A7-A0DC", "TCGA-A7-A13G", "TCGA-CZ-4865", "TCGA-44-2668"], "Cancer type": ["BRCA", "BRCA", "KIRC", "LUAD"],
                    "PanCanAtlas annotation": ["Tumour is DCIS, not invasive (protocol)", "Normal: not normal breast epithelial tissue",
                                               "Prior malignancy treated with taxol and radiation", "Pathology: NSCLC, NOS – to be excluded"],
                    "Decision": ["Kept (author decision 2026-11-09)"] * 4})

tc = rd(R + "telomere_coverage_by_type.tsv").rename(columns={"project": "Cancer type", "tumour_TL": "Tumour TL (n)", "ratio_T/N": "T/N ratio (n)",
     "ratio_diff_centre": "T and N from different centres (n)", "median_log2_TN": "Median log2(T/N)", "frac_T_shorter": "Fraction tumour shorter"})
tlib = tl.dropna(subset=["log2_TN"]).groupby("lib_ratio").log2_TN.agg(n="count", median="median", sd="std",
        frac_T_shorter=lambda s: (s < 0).mean()).reset_index().rename(columns={"lib_ratio": "Sequencing type", "n": "Pairs (n)",
        "median": "Median log2(T/N)", "sd": "SD log2(T/N)", "frac_T_shorter": "Fraction tumour shorter"}).round(2)

P6o = P6.rename(columns={"model": "Model", "clock": "Clock", "n": "n", "beta": "β (proliferation)", "se": "SE", "p": "p", "q": "q (BH)",
                         "beta_dEpi": "β (ΔEpi)", "p_dEpi": "p (ΔEpi)", "r2": "R²"})
W6o = W6.rename(columns={"type": "Cancer type", "clock": "Clock", "n": "n", "beta": "β (proliferation)", "p": "p", "q": "q (BH)"})
Y6o = Y6.rename(columns={"type": "Cancer type", "clock": "Clock", "n": "n", "rho": "Spearman ρ (tumour clock vs yang9)"})
B8o = B8.assign(partial_R2=lambda d: (100 * d.partial_R2).round(2)).rename(columns={"scope": "Model scope", "clock": "Clock", "n": "n",
      "block": "Block", "partial_R2": "Partial R² (%)", "p": "p (F test)", "q": "q (BH)"})
I7o = I7.assign(Status=np.where(I7.outcome.str.startswith(("HNSC", "LIHC")), "Dropped: perfect separation / non-convergence", "Reported")).rename(
      columns={"outcome": "Outcome", "clock": "Clock", "n": "n", "events": "Events (n)", "beta_dclock": "β (Δclock)", "LR": "LR statistic", "p": "p (LR, 1 df)", "q": "q (BH)"})

COV = "Covariates: tumour purity (CPE), ΔEpi, ΔFib, ΔIC (EpiDISH; LIHC composition terms neutralised), signed global Δβ, age, sex, cancer type; continuous variables z-scored within type; HC3 SEs; BH within table."
TABLES = [
 ("S-C8.1", "Table S-C8.1. Data availability per cancer type for the proliferation analyses", "Counts of the 506 audited tumour–normal pairs (13 types) with each data type.", [S1]),
 ("S-C8.2", "Table S-C8.2. Proliferation scores: definitions, agreement and tumour–normal contrast", "Agreement = Spearman ρ within type, tumours only; contrast = fraction of complete pairs with tumour score > normal.", [defs, ag, tn]),
 ("S-C8.3", "Table S-C8.3. Genomic proxies per cancer type and QC-flagged pairs", "ABSOLUTE ploidy and genome doubling (WGD ≥1 doubling); aneuploidy score from Taylor et al. 2018. All ABSOLUTE call types retained.", [gc, qcf]),
 ("S-C8.4", "Table S-C8.4. Telomere length (Barthel et al. 2017): coverage and measurement quality", "Ratio = log2(tumour/solid-normal TL) from the same library type; blood normals not used. Exome-based estimates (WXS) are noisier (limitation).", [tc, tlib]),
 ("S-C8.5", "Table S-C8.5. Pooled models: Δclock ~ tumour proliferation", COV, [P6o]),
 ("S-C8.6", "Table S-C8.6. Within-type models and tumour-only replication of Yang et al. 2016", COV + " Subtype added for BRCA (PAM50) and KIRC.", [W6o, Y6o]),
 ("S-C8.7", "Table S-C8.7. Variance in Δclock explained by each block (drop-one partial R²)", "Blocks: proliferation; composition (CPE, ΔEpi, ΔFib, ΔIC); global Δβ; genomic (aneuploidy, WGD); telomere (ratio + library type). WGS+LPS scope is underpowered and not interpreted.", [B8o]),
 ("S-C8.8", "Table S-C8.8. Information added by Δclock beyond proliferation and global Δβ", "Likelihood-ratio test of adding Δclock to outcome ~ yang9 + CPE + ΔEpi + ΔFib + ΔIC + global Δβ + age + sex (+ type). Ordinal logit (KIRC grade), logit (binary), OLS (aneuploidy).", [I7o]),
]
xl = R + "Supp_Tables_C8.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as w:
    for sheet, title, note, dfs in TABLES:
        row = 3
        pd.DataFrame([[title], [note]]).to_excel(w, sheet_name=sheet, header=False, index=False, startrow=0)
        for d in dfs:
            d.round(4).to_excel(w, sheet_name=sheet, index=False, startrow=row); row += len(d) + 3
        from openpyxl.styles import Font; ws = w.sheets[sheet]; ws["A1"].font = Font(bold=True, size=11)
print(f"supplementary tables -> {xl} ({len(TABLES)} sheets)")
