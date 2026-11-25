#!/usr/bin/env python3
"""C1.1 - mitotic clocks (epiTOC1, epiTOC2, HypoClock, stemTOC) on the rebuilt pairs.

Port of EpiMitClocks (github.com/aet21/EpiMitClocks, R/ *.R); CpG lists/coefficients
exported from its data/dataETOC3.rda and data/epiTOCcpgs3.rda.
  epiTOC1   = mean beta over 385 PCGT CpGs                    (Yang 2016 Genome Biol 17:205)
  epiTOC2   = tnsc = 2*mean((b-b0)/(d*(1-b0))) over 163 CpGs  (Teschendorff 2020 Genome Med 12:56)
              irS = tnsc/age (if age column given)
  HypoClock = 1 - mean beta over 678 solo-WCGW CpGs           (Teschendorff 2020)
  stemTOC   = 95th percentile of beta over 371 CpGs           (Zhu 2024 Nat Commun 15:4211)

Outputs (in --out):
  mitotic_clocks_samples.tsv  one row per pair x tissue (T/N)
  mitotic_clocks_pairs.tsv    one row per pair, T, N and T-N for each clock
  clock_coverage.tsv          CpGs per clock: in clock / on array / non-missing in all samples
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True, help="pairs_rebuilt.tsv")
ap.add_argument("--betadir", required=True, help="folder with beta files")
ap.add_argument("--cpgs", required=True, help="folder with clock_cpgs/*.tsv")
ap.add_argument("--out", default="c1_out")
ap.add_argument("--tcol", default="tumour_file")
ap.add_argument("--ncol", default="normal_file")
ap.add_argument("--idcol", default="patient")
ap.add_argument("--typecol", default="cancer_type")
ap.add_argument("--agecol", default="age")
a = ap.parse_args()

pairs = pd.read_csv(a.pairs, sep="\t")
miss = [c for c in (a.tcol, a.ncol, a.idcol, a.typecol) if c not in pairs.columns]
if miss:
    sys.exit(f"Columns {miss} not in {a.pairs}. Available: {list(pairs.columns)}\n"
             "Pass the right names with --tcol/--ncol/--idcol/--typecol/--agecol.")
has_age = a.agecol in pairs.columns

cp = Path(a.cpgs)
e1 = pd.read_csv(cp / "epiTOC1_cpgs.tsv", sep="\t")["cpg"].tolist()
e2 = pd.read_csv(cp / "epiTOC2_coef.tsv", sep="\t", index_col=0)
hy = pd.read_csv(cp / "HypoClock_cpgs.tsv", sep="\t")["cpg"].tolist()
st = pd.read_csv(cp / "stemTOC_cpgs.tsv", sep="\t")["cpg"].tolist()
clock_cpgs = {"epiTOC1": e1, "epiTOC2": list(e2.index), "HypoClock": hy, "stemTOC": st}
need = set().union(*clock_cpgs.values())

def read_beta(fn):
    p = Path(fn)
    if not p.is_absolute():
        p = Path(a.betadir) / p
    d = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1], names=["cpg", "b"],
                    dtype={"cpg": str}, low_memory=False)
    d = d[d["cpg"].str.startswith(("cg", "ch", "rs"))]          # drops any header line
    on_array = set(d["cpg"])
    d = d[d["cpg"].isin(need)]
    return pd.to_numeric(d.set_index("cpg")["b"], errors="coerce"), on_array

def clocks(b, age):
    r = {}
    v = b.reindex(e1).dropna();             r["epiTOC1"] = v.mean()
    v = b.reindex(e2.index)
    ok = v.notna()
    x = (v[ok] - e2.loc[ok, "beta0"]) / (e2.loc[ok, "delta"] * (1 - e2.loc[ok, "beta0"]))
    r["epiTOC2"] = 2 * x.mean()
    r["epiTOC2_irS"] = r["epiTOC2"] / age if age and age > 0 else np.nan
    v = b.reindex(hy).dropna();             r["HypoClock"] = 1 - v.mean()
    v = b.reindex(st).dropna();             r["stemTOC"] = np.quantile(v, 0.95) if len(v) else np.nan
    return r

rows, nonmiss, array_cpgs = [], {k: None for k in clock_cpgs}, None
for i, pr in pairs.iterrows():
    age = pr[a.agecol] if has_age else None
    for tis, col in (("T", a.tcol), ("N", a.ncol)):
        b, on_arr = read_beta(pr[col])
        array_cpgs = on_arr if array_cpgs is None else array_cpgs & on_arr
        for k, cs in clock_cpgs.items():
            present = set(b.reindex(cs).dropna().index)
            nonmiss[k] = present if nonmiss[k] is None else nonmiss[k] & present
        rows.append({"patient": pr[a.idcol], "cancer_type": pr[a.typecol], "tissue": tis,
                     "age": age, "file": pr[col], **clocks(b, age)})
    if (i + 1) % 50 == 0:
        print(f"{i+1}/{len(pairs)} pairs", file=sys.stderr)

out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
s = pd.DataFrame(rows)
s.to_csv(out / "mitotic_clocks_samples.tsv", sep="\t", index=False)

vals = ["epiTOC1", "epiTOC2", "epiTOC2_irS", "HypoClock", "stemTOC"]
w = s.pivot(index=["patient", "cancer_type"], columns="tissue", values=vals)
w.columns = [f"{c}_{t}" for c, t in w.columns]
w = w.join(s.groupby(["patient", "cancer_type"])["age"].first())
for c in vals:
    w[f"{c}_delta"] = w[f"{c}_T"] - w[f"{c}_N"]
w.reset_index().to_csv(out / "mitotic_clocks_pairs.tsv", sep="\t", index=False)

cov = pd.DataFrame([{"clock": k, "n_clock": len(cs),
                     "n_on_array": len(set(cs) & array_cpgs),
                     "n_nonmissing_all_samples": len(nonmiss[k])} for k, cs in clock_cpgs.items()])
cov.to_csv(out / "clock_coverage.tsv", sep="\t", index=False)
print(cov.to_string(index=False))
print(f"\n{len(s)} samples, {len(w)} pairs -> {out}/")
