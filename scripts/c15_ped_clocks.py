#!/usr/bin/env python3
"""C1.5 - mitotic clocks (epiTOC1, epiTOC2, HypoClock, stemTOC) per pediatric sample (no pairing needed).
Same clock code as c1_mitotic_clocks_v2.py. One row per row of --map; adds per-sample CpG counts."""
import argparse, sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True, help="Sample_mapping_*.tsv (needs file_path, age_at_diagnosis_years)")
ap.add_argument("--base", default=".", help="folder that file_path is relative to")
ap.add_argument("--cpgs", required=True, help="folder with clock_cpgs/*.tsv")
ap.add_argument("--agecol", default="age_at_diagnosis_years")
ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--out", default="c15_out")
a = ap.parse_args()

m = pd.read_csv(a.map, sep="\t")
cp = Path(a.cpgs)
e1 = pd.read_csv(cp/"epiTOC1_cpgs.tsv", sep="\t")["cpg"].tolist()
e2 = pd.read_csv(cp/"epiTOC2_coef.tsv", sep="\t", index_col=0)
hy = pd.read_csv(cp/"HypoClock_cpgs.tsv", sep="\t")["cpg"].tolist()
st = pd.read_csv(cp/"stemTOC_cpgs.tsv", sep="\t")["cpg"].tolist()
CL = {"epiTOC1": e1, "epiTOC2": list(e2.index), "HypoClock": hy, "stemTOC": st}
need = set().union(*CL.values())

def one(args):
    fp, age = args
    try:
        d = pd.read_csv(Path(a.base)/fp, sep="\t", header=None, usecols=[0, 1], names=["cpg", "b"],
                        dtype={"cpg": str}, low_memory=False)
        d = d[d["cpg"].str.startswith(("cg", "ch", "rs"))]
        n_probes = len(d)
        b = pd.to_numeric(d[d["cpg"].isin(need)].drop_duplicates("cpg").set_index("cpg")["b"], errors="coerce")
        r = {"n_probes": n_probes}
        r["epiTOC1"] = b.reindex(e1).dropna().mean()
        v = b.reindex(e2.index); ok = v.notna()
        r["epiTOC2"] = 2*((v[ok]-e2.loc[ok, "beta0"])/(e2.loc[ok, "delta"]*(1-e2.loc[ok, "beta0"]))).mean()
        r["epiTOC2_irS"] = r["epiTOC2"]/age if age and age > 0 else np.nan
        r["HypoClock"] = 1 - b.reindex(hy).dropna().mean()
        v = b.reindex(st).dropna(); r["stemTOC"] = np.quantile(v, .95) if len(v) else np.nan
        for k, cs in CL.items():
            r[f"nCpG_{k}"] = int(b.reindex(cs).notna().sum())
        return r
    except Exception as e:
        print("FAILED", fp, e, file=sys.stderr); return {}

if __name__ == "__main__":
    jobs = list(zip(m["file_path"], pd.to_numeric(m[a.agecol], errors="coerce")))
    with Pool(a.jobs) as p:
        res = []
        for i, r in enumerate(p.imap(one, jobs, chunksize=4), 1):
            res.append(r)
            if i % 200 == 0: print(f"{i}/{len(jobs)}", file=sys.stderr)
    out = Path(a.out); out.mkdir(exist_ok=True, parents=True)
    s = pd.concat([m.reset_index(drop=True), pd.DataFrame(res)], axis=1)
    s.to_csv(out/"ped_mitotic_clocks_samples.tsv", sep="\t", index=False)
    cov = pd.DataFrame([{"clock": k, "n_clock": len(c), "per_sample_min": s[f"nCpG_{k}"].min(),
                         "per_sample_median": s[f"nCpG_{k}"].median(),
                         "n_lt80pct": int((s[f"nCpG_{k}"] < .8*len(c)).sum())} for k, c in CL.items()])
    cov.to_csv(out/"ped_clock_coverage.tsv", sep="\t", index=False)
    print(cov.to_string(index=False)); print("n_probes per file: min/median/max",
          s.n_probes.min(), s.n_probes.median(), s.n_probes.max()); print(f"{len(s)} samples -> {out}/")
