#!/usr/bin/env python3
"""Dump betas at all clock CpGs (union of 4 clocks) per pediatric sample -> wide matrix, so clocks can be
recomputed on any common CpG set (arrays differ: 450K / EPIC / EPICv2)."""
import argparse, sys
from pathlib import Path
from multiprocessing import Pool
import pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True); ap.add_argument("--base", default=".")
ap.add_argument("--cpgs", required=True); ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--out", default="c15_out/ped_clock_cpg_betas.tsv.gz")
a = ap.parse_args()
cp = Path(a.cpgs)
need = sorted(set(pd.read_csv(cp/"epiTOC1_cpgs.tsv", sep="\t")["cpg"]) | set(pd.read_csv(cp/"epiTOC2_coef.tsv", sep="\t", index_col=0).index)
              | set(pd.read_csv(cp/"HypoClock_cpgs.tsv", sep="\t")["cpg"]) | set(pd.read_csv(cp/"stemTOC_cpgs.tsv", sep="\t")["cpg"]))
def one(fp):
    try:
        d = pd.read_csv(Path(a.base)/fp, sep="\t", header=None, usecols=[0, 1], names=["cpg", "b"], dtype={"cpg": str}, low_memory=False)
        d = d[d["cpg"].isin(need)].drop_duplicates("cpg").set_index("cpg")["b"]
        return pd.to_numeric(d, errors="coerce").reindex(need)
    except Exception as e:
        print("FAILED", fp, e, file=sys.stderr); return pd.Series(index=need, dtype=float)
if __name__ == "__main__":
    m = pd.read_csv(a.map, sep="\t").dropna(subset=["file_path"])
    with Pool(a.jobs) as p: res = p.map(one, m["file_path"], chunksize=4)
    M = pd.DataFrame(res, index=m["file_uuid"].values).round(4)
    Path(a.out).parent.mkdir(exist_ok=True, parents=True); M.to_csv(a.out, sep="\t", compression="gzip")
    print(M.shape, "->", a.out)
