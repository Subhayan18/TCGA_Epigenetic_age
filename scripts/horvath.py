import pandas as pd, numpy as np
D = "/home/chattopa/data_storage/TCGA_GDC_rebuilt"
w = pd.read_csv("Horvath1.csv", index_col=0).iloc[:, 0]
g = pd.read_csv("gdc_barcodes.tsv", sep="\t"); g.columns = ["aliquot", "project", "pid", "file", "uuid"]
u = dict(zip(g.file, g.uuid)); p = pd.read_csv("pairs_rebuilt.tsv", sep="\t")
rd = lambda i: pd.read_csv(f"{D}/{i}.txt", sep="\t", header=None, index_col=0).iloc[:, 0].reindex(w.index)
M = pd.DataFrame({(r.pid, s): rd(u[getattr(r, s + "_file")]) for r in p.itertuples() for s in "tn"})
nm = M.isna().sum(); print("CpGs absent from all files:", int(M.isna().all(axis=1).sum()), "| max missing/sample:", int(nm.max()))
M = M.where(M.notna(), M.mean(axis=1), axis=0)
x = M.mul(w, axis=0).sum() + 0.695507258
a = pd.Series(np.where(x < 0, 21 * np.exp(x) - 1, 21 * x + 20), index=M.columns)
o = pd.DataFrame({"patient_id": p.type + "_" + p.pid.str[5:], "DNAmAge_Normal": [a[(i, "n")] for i in p.pid], "DNAmAge_Tumour": [a[(i, "t")] for i in p.pid],
                  "Chronological_Age": p.age, "n_missing_N": [nm[(i, "n")] for i in p.pid], "n_missing_T": [nm[(i, "t")] for i in p.pid]})
o.to_csv("paired_results_rebuilt.csv", index=False); print(o.shape); print(o.describe().round(1))
