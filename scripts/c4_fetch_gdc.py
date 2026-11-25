#!/usr/bin/env python3
"""C4.1/C4.3 inputs from GDC: sex per case + ASCAT allele-specific CN segments
(primary tumour, open access), then a per-patient X copy-number summary.
Run from Meth_TCGA_GSE213478/:  python3 scripts/c4_fetch_gdc.py"""
import json, os, re, time, urllib.request
import numpy as np, pandas as pd

API = "https://api.gdc.cancer.gov"
RAW = "/home/chattopa/data_storage/TCGA_GDC_ASCAT"          # raw segment files, next to the betas
PAIRS = next(p for p in ["data/pairs/pairs_gdc.tsv", "pairs_gdc.tsv"] if os.path.exists(p))
for d in [RAW, "data/gdc_clinical", "data/gdc_ascat", "results/c4_xchrom"]:
    os.makedirs(d, exist_ok=True)

def post(ep, payload):
    req = urllib.request.Request(f"{API}/{ep}", json.dumps(payload).encode(), {"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))["data"]["hits"]

def get(fid, dest, tries=3):
    for i in range(tries):
        try:
            urllib.request.urlretrieve(f"{API}/data/{fid}", dest + ".part"); os.replace(dest + ".part", dest); return
        except Exception as e:
            print("retry", fid, e); time.sleep(5 * (i + 1))
    raise RuntimeError(f"download failed: {fid}")

ids = sorted(set(re.findall(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", open(PAIRS).read())))
print(len(ids), "cases from", PAIRS)
f_in = lambda field, v: {"op": "in", "content": {"field": field, "value": v}}

# 1) clinical: sex, race, age
hits = post("cases", {"filters": {"op": "and", "content": [f_in("submitter_id", ids), f_in("project.program.name", ["TCGA"])]},
                      "fields": "submitter_id,project.project_id,demographic.gender,demographic.sex_at_birth,demographic.race,diagnoses.age_at_diagnosis",
                      "size": 5000})
clin = pd.DataFrame([{"pid": h["submitter_id"], "project": h["project"]["project_id"],
                      "sex": (h.get("demographic") or {}).get("sex_at_birth") or (h.get("demographic") or {}).get("gender"), "race": (h.get("demographic") or {}).get("race"),
                      "age_dx_days": ((h.get("diagnoses") or [{}])[0]).get("age_at_diagnosis")} for h in hits])
clin.to_csv("data/gdc_clinical/case_clinical.tsv", sep="\t", index=False)
print("clinical:", len(clin), "cases;", clin["sex"].value_counts(dropna=False).to_dict())

# 2) ASCAT allele-specific segment files, one per case (primary tumour; 01A first, ASCAT2 before ASCAT3)
hits = post("files", {"filters": {"op": "and", "content": [f_in("cases.submitter_id", ids),
                      f_in("data_type", ["Allele-specific Copy Number Segment"]), f_in("access", ["open"])]},
                      "fields": "file_id,file_name,analysis.workflow_type,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type",
                      "size": 10000})
rows = []
for h in hits:
    for c in h["cases"]:
        for s in c.get("samples", []):
            if s.get("sample_type") == "Primary Tumor":
                rows.append({"pid": c["submitter_id"], "ascat_sample": s["submitter_id"], "file_id": h["file_id"],
                             "file_name": h["file_name"], "workflow": h.get("analysis", {}).get("workflow_type")})
man = pd.DataFrame(rows)
print("ASCAT files (primary tumour):", len(man), "| workflows:", man["workflow"].value_counts().to_dict())
man["r"] = (~man.ascat_sample.str.endswith("01A")).astype(int) * 10 + man.workflow.map({"ASCAT2": 0, "ASCAT3": 1}).fillna(2)
man = man.sort_values(["pid", "r", "file_name"]).drop_duplicates("pid").drop(columns="r")
man.to_csv("data/gdc_ascat/ascat_manifest.tsv", sep="\t", index=False)
print("cases with ASCAT:", len(man), "of", len(ids))

# 3) download + per-patient X summary (whole-chromosome; build-independent)
AUTO = [str(i) for i in range(1, 23)]
out = []
for _, m in man.iterrows():
    dest = f"{RAW}/{m.file_id}.tsv"
    if not os.path.exists(dest):
        get(m.file_id, dest)
    s = pd.read_csv(dest, sep="\t")
    s["chrom"] = s["Chromosome"].astype(str).str.replace("chr", "", regex=False)
    s["len"] = s["End"] - s["Start"] + 1
    a, x = s[s.chrom.isin(AUTO)], s[s.chrom == "X"]
    wm = lambda d, v: np.average(d[v], weights=d["len"]) if len(d) else np.nan
    loh = lambda d: d.loc[(d.Minor_Copy_Number == 0) & (d.Copy_Number > 0), "len"].sum() / d["len"].sum() if len(d) else np.nan
    out.append({"pid": m.pid, "ascat_sample": m.ascat_sample, "workflow": m.workflow,
                "auto_mean_cn": wm(a, "Copy_Number"), "x_mean_cn": wm(x, "Copy_Number"),
                "x_loh_frac": loh(x), "auto_loh_frac": loh(a)})
x = pd.DataFrame(out)
x["x_rel_cn"] = x.x_mean_cn / x.auto_mean_cn               # ~1 female diploid X, ~0.5 male
res = clin.merge(x, on="pid", how="left")
res.to_csv("results/c4_xchrom/x_cna_by_patient.tsv", sep="\t", index=False)
print(res.groupby("sex")[["x_rel_cn", "x_loh_frac"]].describe().round(2).T.to_string())
