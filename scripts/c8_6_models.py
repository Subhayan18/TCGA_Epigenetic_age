#!/usr/bin/env python3
"""C8.6-C8.9 models: what does the tumour-normal mitotic-clock advance track?
Run from Meth_TCGA_GSE213478/:
    python3 scripts/c8_6_models.py --check     # only report which input columns were matched
    python3 scripts/c8_6_models.py             # all models
Needs pandas, numpy, scipy, statsmodels. Outputs -> results/c8_prolif/models/
Design (agreed 2026-11-09):
  C8.6 dclock ~ tumour proliferation (yang9) + CPE + dEpi + dFib + dIC + signed mean dbeta + age + sex [+ type FE];
       pooled + within type (+PAM50 BRCA, mRNA subtype KIRC); sensitivity dRNA (complete T+N pairs), G2M;
       supp: tumour-only clock ~ yang9 (Yang 2016 replication)
  C8.7 incremental value of dclock beyond proliferation: KIRC grade (ordinal), HNSC/LIHC G3 vs G1-2,
       late stage pooled KIRC/BRCA/KIRP/THCA, WGD pooled, aneuploidy pooled (excl THCA, PRAD); LR tests
  C8.8 block contributions (drop-one partial R2): proliferation | composition | global dbeta | genomic | telomere
  C8.9 sensitivity: CPE>=0.6, without global dbeta, without low-coverage pairs, telomere WGS+LPS only
All continuous variables z-scored within cancer type; HC3 SEs; BH within each table.
"""
import argparse, os, re, sys, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr, chi2, f as fdist

ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); A = ap.parse_args()
O = "results/c8_prolif/models"; os.makedirs(O, exist_ok=True)
log = lambda *a: print(*a, flush=True)
TYPES13 = ["KIRC","BRCA","PRAD","HNSC","THCA","LIHC","LUSC","KIRP","UCEC","BLCA","LUAD","COAD","ESCA"]
CLOCKS = ["epiTOC1", "epiTOC2", "HypoClock", "stemTOC", "Horvath"]
LOWCOV = ["TCGA-44-6145", "TCGA-B0-5402"]

def read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t", low_memory=False)
def pid_index(df, name):
    for c in df.columns:
        s = df[c].astype(str)
        if s.str.match(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}").mean() > 0.9:
            df = df.copy(); df.index = s.str[:12]; return df[~df.index.duplicated()]
    for c in df.columns:  # e.g. BLCA_BL-A13J -> TCGA-BL-A13J
        s = df[c].astype(str)
        if s.str.contains(r"[A-Z0-9]{2}-[A-Z0-9]{4}").mean() > 0.9:
            df = df.copy(); df.index = "TCGA-" + s.str.extract(r"([A-Z0-9]{2}-[A-Z0-9]{4})")[0]; return df[~df.index.duplicated()]
    sys.exit(f"[{name}] no TCGA barcode column. Header: {list(df.columns)}")
def find(df, pat, name, required=True):
    hits = [c for c in df.columns if re.search(pat, str(c), re.I)]
    if hits: return hits[0]
    if required: sys.exit(f"[{name}] no column matching /{pat}/. Header: {list(df.columns)}")
    return None

# ---------------- assemble patient table ----------------
M = pid_index(read("c1_out/mitotic_clocks_pairs.tsv"), "mitotic")
H = pid_index(read("pair_table_corrected.csv"), "horvath")
F = pid_index(read("c2_out/fractions_pairs.tsv"), "fractions")
G = pid_index(read("results/c3_phylo/hetero_v2/hetero_clusters_patients.tsv"), "hetero")
CC = pid_index(read("data/gdc_clinical/case_clinical.tsv"), "clinical")
R = "results/c8_prolif/"
PR = read(R + "prolif_scores_patients.tsv").set_index("patient")
GE = read(R + "genomic_patients.tsv").set_index("patient")
TL = read(R + "telomere_patients.tsv").set_index("patient")
CL = read(R + "clinical_patients.tsv").set_index("patient")

D = pd.DataFrame(index=PR.index)
used = {}
for k in CLOCKS[:4]:
    t, n = find(M, rf"^{k}.*_T$", "mitotic"), find(M, rf"^{k}.*_N$", "mitotic")
    D[f"{k}_T"], D[f"d_{k}"] = M[t], M[t] - M[n]; used[k] = (t, n)
t, n = find(H, r"DNAmAge_Tumou?r", "horvath"), find(H, r"DNAmAge_Normal", "horvath")
D["Horvath_T"], D["d_Horvath"] = H[t], H[t] - H[n]; used["Horvath"] = (t, n)
for c in ("Epi", "Fib", "IC"):
    t = find(F, rf"^{c}(_|\.)?T$|^T(_|\.){c}$|{c}_tum", "fractions")
    n = find(F, rf"^{c}(_|\.)?N$|^N(_|\.){c}$|{c}_norm", "fractions")
    D[f"d{c}"] = F[t] - F[n]; used[c] = (t, n)
cpe = find(F, r"^CPE$|cpe", "fractions", required=False)
if cpe: D["CPE"] = F[cpe]; used["CPE"] = cpe
else:
    P = read("data/reference/c2_ref/tcga_purity_aran2015.tsv"); P = pid_index(P, "purity")
    D["CPE"] = pd.to_numeric(P[find(P, r"^CPE$", "purity")].astype(str).str.replace(",", "."), errors="coerce"); used["CPE"] = "aran2015"
gb = find(G, r"signed.*mean|mean.*d(elta)?_?beta|^mean_d", "hetero"); D["gdb"] = G[gb]; used["global_dbeta"] = gb
D["age"] = pd.to_numeric(CC[find(CC, r"age", "clinical")], errors="coerce")
D["male"] = (CC[find(CC, r"sex|gender", "clinical")].astype(str).str.lower() == "male").astype(float)
D["type"] = PR.project
D = D.join(PR[["yang9_T", "g2m_T", "yang9_D", "g2m_D"]])
D = D.join(GE[["ploidy", "wgd", "aneuploidy"]])
D = D.join(TL[["log2_TN", "lib_ratio"]])
D = D.join(CL[["grade", "stage_late", "subtype_mRNA"]])
D = D[D.type.isin(TYPES13)]
log("matched columns:"); [log(f"  {k}: {v}") for k, v in used.items()]
log(f"patients: {len(D)}; non-missing: " + ", ".join(f"{c} {D[c].notna().sum()}" for c in
    ["d_stemTOC", "d_Horvath", "dEpi", "CPE", "gdb", "age", "yang9_T", "yang9_D", "aneuploidy", "wgd", "log2_TN", "grade"]))
if A.check: sys.exit(0)

# ---------------- preprocessing ----------------
w = D.log2_TN.quantile([.01, .99]); D["tl"] = D.log2_TN.clip(*w)
D.loc[D.type == "LIHC", ["dEpi", "dFib", "dIC"]] = 0.0     # LIHC deconvolution invalid (C2.1): composition terms neutralised
Z = D.copy()
cont = [f"d_{k}" for k in CLOCKS] + [f"{k}_T" for k in CLOCKS] + ["yang9_T", "g2m_T", "yang9_D", "g2m_D", "CPE",
        "dEpi", "dFib", "dIC", "gdb", "age", "ploidy", "aneuploidy", "tl"]
for c in cont:
    Z[c] = D.groupby("type")[c].transform(lambda v: (v - v.mean()) / v.std(ddof=0) if v.std(ddof=0) > 0 else v * 0)
BASE = "CPE + dEpi + dFib + dIC + gdb + age + male"

def ols(f, d):
    try: return smf.ols(f, d).fit(cov_type="HC3")
    except Exception as e: print("  ols failed:", f, "|", e); return None
def bh(df, col="p"):
    m = df[col].notna(); df.loc[m, "q"] = multipletests(df.loc[m, col], method="fdr_bh")[1]; return df

# ---------------- C8.6 ----------------
rows = []
def c86(label, d, prolif, base=BASE, fe=True, types=None):
    for k in CLOCKS:
        f = f"d_{k} ~ {prolif} + {base}" + (" + C(type)" if fe else "")
        x = d.dropna(subset=[f"d_{k}", prolif] + [v.strip() for v in base.split("+") if v.strip()])
        if types: x = x[x.type.isin(types)]
        r = ols(f, x)
        if r is not None and prolif in r.params:
            rows.append(dict(model=label, clock=k, n=int(r.nobs), beta=r.params[prolif], se=r.bse[prolif], p=r.pvalues[prolif],
                             beta_dEpi=r.params.get("dEpi", np.nan), p_dEpi=r.pvalues.get("dEpi", np.nan), r2=r.rsquared))
c86("pooled_main_yang9", Z, "yang9_T")
c86("pooled_G2M", Z, "g2m_T")
c86("pooled_dRNA_yang9", Z, "yang9_D")
c86("sens_CPE>=0.6", Z[D.CPE >= 0.6], "yang9_T")
c86("sens_no_global_dbeta", Z, "yang9_T", base=BASE.replace("gdb + ", ""))
c86("sens_drop_lowcov", Z.drop(LOWCOV, errors="ignore"), "yang9_T")
c86("sens_noCPE_incl_ESCA", Z, "yang9_T", base=BASE.replace("CPE + ", ""))
c86("sens_excl_KIRC", Z[Z.type != "KIRC"], "yang9_T")
T6 = bh(pd.DataFrame(rows)); T6.to_csv(f"{O}/C8_6_pooled.tsv", sep="\t", index=False)
log("\nC8.6 pooled: dclock ~ proliferation (+ composition, global dbeta, age, sex, type)\n" +
    T6.pivot(index="model", columns="clock", values="beta").round(2).to_string())
log("q-values:\n" + T6.pivot(index="model", columns="clock", values="q").map(lambda v: f"{v:.1e}").to_string())

rows = []
for t, d in Z.groupby("type"):
    base = BASE + (" + C(subtype_mRNA)" if t in ("BRCA", "KIRC") else "")
    if t == "ESCA": base = base.replace("CPE + ", "")
    for k in CLOCKS:
        x = d.dropna(subset=[f"d_{k}", "yang9_T", "age"])
        if t in ("BRCA", "KIRC"): x = x.dropna(subset=["subtype_mRNA"])
        if len(x) < 12: continue
        r = ols(f"d_{k} ~ yang9_T + {base}", x)
        if r is not None:
            rows.append(dict(type=t, clock=k, n=int(r.nobs), beta=r.params["yang9_T"], p=r.pvalues["yang9_T"]))
TW = bh(pd.DataFrame(rows)); TW.to_csv(f"{O}/C8_6_within_type.tsv", sep="\t", index=False)
log("\nC8.6 within type, beta for yang9 (* q<0.05):\n" + TW.assign(
    b=TW.beta.round(2).astype(str) + np.where(TW.q < 0.05, "*", "")).pivot(index="type", columns="clock", values="b").to_string())

rows = [dict(type=t, clock=k, n=d[[f"{k}_T", "yang9_T"]].dropna().shape[0],
             rho=spearmanr(d[f"{k}_T"], d.yang9_T, nan_policy="omit")[0]) for t, d in D.groupby("type") for k in CLOCKS]
S = pd.DataFrame(rows); S.to_csv(f"{O}/C8_6_supp_tumour_only_yang.tsv", sep="\t", index=False)
log("\nsupp (Yang 2016 replication), tumour clock vs yang9, within-type Spearman:\n" +
    S.pivot(index="type", columns="clock", values="rho").round(2).to_string())

# ---------------- C8.8 block contributions ----------------
BLOCKS = {"proliferation": ["yang9_T"], "composition": ["CPE", "dEpi", "dFib", "dIC"], "global_dbeta": ["gdb"],
          "genomic": ["aneuploidy", "wgd"], "telomere": ["tl", "C(lib_ratio)"]}
rows = []
for scope, blocks in [("no_telomere", [b for b in BLOCKS if b != "telomere"]), ("with_telomere", list(BLOCKS)),
                      ("telomere_WGS_LPS", list(BLOCKS))]:
    x0 = Z if scope != "telomere_WGS_LPS" else Z[Z.lib_ratio.isin(["WGS", "LPS"])]
    for k in CLOCKS:
        terms = sum((BLOCKS[b] for b in blocks), []) + ["age", "male", "C(type)"]
        if scope == "telomere_WGS_LPS": terms = [t for t in terms if t != "C(lib_ratio)"]
        need = [t for t in terms if not t.startswith("C(")] + [f"d_{k}"] + (["lib_ratio"] if "C(lib_ratio)" in terms else [])
        x = x0.dropna(subset=need)
        if len(x) < 30: continue
        for b in blocks:
            red_terms = [t for t in terms if t not in BLOCKS[b]]
            red = smf.ols(f"d_{k} ~ " + " + ".join(red_terms), x).fit()
            fu = smf.ols(f"d_{k} ~ " + " + ".join(terms), x).fit()
            df_num = fu.df_model - red.df_model
            F_ = ((red.ssr - fu.ssr) / df_num) / (fu.ssr / fu.df_resid) if df_num > 0 else np.nan
            rows.append(dict(scope=scope, clock=k, n=int(fu.nobs), block=b, partial_R2=(red.ssr - fu.ssr) / red.ssr,
                             p=fdist.sf(F_, df_num, fu.df_resid) if df_num > 0 else np.nan))
B8 = pd.DataFrame(rows)
B8 = pd.concat([bh(g.copy()) for _, g in B8.groupby("scope")]); B8.to_csv(f"{O}/C8_8_block_partialR2.tsv", sep="\t", index=False)
for s, g in B8.groupby("scope", sort=False):
    log(f"\nC8.8 partial R2 by block [{s}, n={g.n.max()}] (* q<0.05):\n" + g.assign(
        v=(100 * g.partial_R2).round(1).astype(str) + "%" + np.where(g.q < 0.05, "*", "")).pivot(index="block", columns="clock", values="v").to_string())

# ---------------- C8.7 incremental value ----------------
def lr(f0, f1, d, kind):
    try:
        if kind == "logit":
            m0, m1 = smf.logit(f0, d).fit(disp=0), smf.logit(f1, d).fit(disp=0)
        elif kind == "ordinal":
            from statsmodels.miscmodels.ordinal_model import OrderedModel
            m0 = OrderedModel.from_formula(f0.replace(" + 1", ""), d, distr="logit").fit(method="bfgs", disp=0)
            m1 = OrderedModel.from_formula(f1.replace(" + 1", ""), d, distr="logit").fit(method="bfgs", disp=0)
        else:
            m0, m1 = smf.ols(f0, d).fit(), smf.ols(f1, d).fit()
        v = f1.split("+")[-1].strip(); stat = 2 * (m1.llf - m0.llf); return stat, chi2.sf(stat, 1), m1.params.get(v)
    except Exception as e:
        return np.nan, np.nan, np.nan
rows = []
Z["G3"] = Z.grade.map({"G1": 0, "G2": 0, "G3": 1, "G4": 1})
Z["grade_ord"] = pd.Categorical(Z.grade, categories=["G2", "G3", "G4"], ordered=True)
specs = [("KIRC grade (ordinal)", Z[Z.type == "KIRC"], "grade_ord", "ordinal", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male"),
         ("HNSC G3 vs G1-2", Z[Z.type == "HNSC"], "G3", "logit", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male"),
         ("LIHC G3 vs G1-2", Z[Z.type == "LIHC"], "G3", "logit", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male"),
         ("late stage pooled", Z[Z.type.isin(["KIRC", "BRCA", "KIRP", "THCA"])], "stage_late", "logit", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male + C(type)"),
         ("late stage pooled excl KIRC", Z[Z.type.isin(["BRCA", "KIRP", "THCA"])], "stage_late", "logit", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male + C(type)"),
         ("WGD pooled", Z, "wgd", "logit", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male + C(type)"),
         ("aneuploidy pooled (excl THCA, PRAD)", Z[~Z.type.isin(["THCA", "PRAD"])], "aneuploidy", "ols", "yang9_T + CPE + dEpi + dFib + dIC + gdb + age + male + C(type)")]
for name, d, y, kind, base in specs:
    for k in CLOCKS:
        x = d.dropna(subset=[y, f"d_{k}", "yang9_T", "CPE", "age"])
        f0 = f"{y} ~ {base}"; f1 = f"{f0} + d_{k}"
        stat, p, b = lr(f0, f1, x, kind)
        rows.append(dict(outcome=name, clock=k, n=len(x), events=(x[y].sum() if kind == "logit" else np.nan),
                         beta_dclock=b, LR=stat, p=p))
T7 = pd.DataFrame(rows); T7 = pd.concat([bh(g.copy()) for _, g in T7.groupby("outcome")])
T7.to_csv(f"{O}/C8_7_incremental.tsv", sep="\t", index=False)
log("\nC8.7 incremental value of dclock beyond proliferation (LR test p; * q<0.05):\n" + T7.assign(
    v=T7.p.map(lambda v: f"{v:.3f}") + np.where(T7.q < 0.05, "*", "")).pivot(index="outcome", columns="clock", values="v").to_string())
log(f"\noutputs in {O}/")
