#!/usr/bin/env python3
"""C8 data download (C8.1-C8.5). Run from Meth_TCGA_GSE213478/:
    python3 scripts/c8_1_download.py                # everything
    python3 scripts/c8_1_download.py --only rna     # one block: rna | tables | genesets | clinical
Resumable: existing files with correct size/md5 are skipped. Progress printed to terminal.
Outputs: data/c8/ (tables, manifests, matrix); raw RNA files in --rna_dir.
"""
import argparse, csv, gzip, hashlib, json, os, re, sys, time, urllib.request, urllib.error

API = "https://api.gdc.cancer.gov"
TYPES13 = ["KIRC","BRCA","PRAD","HNSC","THCA","LIHC","LUSC","KIRP","UCEC","BLCA","LUAD","COAD","ESCA"]

# External tables (UUIDs from GDC publication pages: pancanatlas, panimmune, pancan-aneuploidy)
GDC_TABLES = {
    "absolute_purity_ploidy.tsv":   "4f277128-f793-4354-a13d-30cc7fe9f6b5",  # TCGA_mastercalls.abs_tables_JSedit.fixed.txt (ploidy, Genome doublings)
    "aneuploidy_arm_calls.tsv":     "4c35f34f-b0f3-4891-8794-4840dd748aad",  # PANCAN_ArmCallsAndAneuploidyScore_092817.txt (Taylor 2018)
    "tcga_cdr_clinical.xlsx":       "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",  # TCGA-CDR (Liu 2018): grade, stage
    "tcga_subtypes.tsv":            "0f31b768-7f67-4fc4-abc3-06ac5bd90bf0",  # TCGASubtype.20170308.tsv (PAM50 etc.)
    "panimmune_scores160.tsv.gz":   "80a82092-161d-4615-9d96-e858f113618d",  # Scores_160_Signatures (Thorsson 2018; proliferation module)
    "panimmune_geneset_defs.xlsx":  "9b174979-fe97-48bc-9e97-9384b0519f03",
    "sample_quality_annotations.tsv":"1a7d7be8-675d-4e60-a105-19d4121bdebf", # PanCanAtlas sample QC/exclusion flags
}
OTHER_URLS = {  # not on GDC; if a download fails, get it by hand (URL printed) into data/c8/
    "barthel2017_telomere_ST1.xlsx":
        "https://static-content.springer.com/esm/art%3A10.1038%2Fng.3781/MediaObjects/41588_2017_BFng3781_MOESM2_ESM.xlsx",
}
GENESETS = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"]
MSIGDB = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName={}&fileType=json"


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def http(url, data=None, tries=5):
    hdr = {"Content-Type": "application/json"} if data else {}
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdr)
            return urllib.request.urlopen(req, timeout=120)
        except (urllib.error.URLError, TimeoutError) as e:
            log(f"  retry {i+1}/{tries}: {e}"); time.sleep(5 * (i + 1))
    raise SystemExit(f"FAILED: {url}")


def fetch(url, out, md5=None, size=None, tries=3):
    if os.path.exists(out) and (size is None or os.path.getsize(out) == size):
        if md5 is None or hashlib.md5(open(out, "rb").read()).hexdigest() == md5:
            return "skip"
    tmp = out + ".part"
    for i in range(tries):
        with http(url) as r, open(tmp, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b: break
                f.write(b)
        got = hashlib.md5(open(tmp, "rb").read()).hexdigest() if md5 else None
        if md5 is None or got == md5:
            os.replace(tmp, out); return "ok"
        log(f"  md5 mismatch {os.path.basename(out)} (size {os.path.getsize(tmp)} vs {size}), retry {i+1}/{tries}")
        time.sleep(5)
    os.remove(tmp); raise SystemExit(f"md5 mismatch: {out}")


def read_pairs(path):
    """Extract patient + methylation tumour/normal sample barcodes from any columns of pairs_gdc.tsv."""
    bc = re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-\d{2}[A-Z]?)?")
    pairs = {}
    with open(path) as f:
        for row in f:
            hits = bc.findall(row)
            if not hits: continue
            pid = hits[0][:12]
            t = next((h[:16] for h in hits if h[13:15] == "01"), None)
            n = next((h[:16] for h in hits if h[13:15] == "11"), None)
            old = pairs.get(pid, (None, None))
            pairs[pid] = (t or old[0], n or old[1])
    return pairs


def gdc_query(endpoint, filters, fields, size=10000):
    body = json.dumps({"filters": filters, "fields": ",".join(fields), "format": "JSON", "size": size}).encode()
    return json.load(http(f"{API}/{endpoint}", data=body))["data"]["hits"]


def block_rna(pairs, outdir, rna_dir):
    log(f"RNA: querying GDC STAR-Counts for {len(pairs)} patients")
    filt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.submitter_id", "value": sorted(pairs)}},
        {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor", "Solid Tissue Normal"]}},
        {"op": "=", "content": {"field": "access", "value": "open"}}]}
    hits = gdc_query("files", filt, ["file_id", "file_name", "md5sum", "file_size", "cases.submitter_id",
                     "cases.project.project_id", "cases.samples.submitter_id", "cases.samples.sample_type",
                     "cases.samples.portions.analytes.aliquots.submitter_id"])
    log(f"  {len(hits)} candidate files")
    cand = {}  # (pid, T/N) -> list
    for h in hits:
        c = h["cases"][0]; s = c["samples"][0]
        side = "T" if s["sample_type"] == "Primary Tumor" else "N"
        aliq = s["portions"][0]["analytes"][0]["aliquots"][0]["submitter_id"]
        cand.setdefault((c["submitter_id"], side), []).append(dict(
            file_id=h["file_id"], md5=h["md5sum"], size=h["file_size"], sample=s["submitter_id"][:16],
            aliquot=aliq, project=c["project"]["project_id"].replace("TCGA-", "")))
    # selection: same sample barcode as methylation (e.g. -01A / -11A) if available, else lowest aliquot barcode
    man = []
    for (pid, side), lst in cand.items():
        meth = pairs[pid][0 if side == "T" else 1]
        same = [x for x in lst if x["sample"] == meth]
        x = sorted(same or lst, key=lambda d: d["aliquot"])[0]
        man.append(dict(patient=pid, side=side, meth_sample=meth, rna_sample_match=bool(same),
                        n_candidates=len(lst), **x))
    man = [m for m in man if m["project"] in TYPES13]
    with open(f"{outdir}/rna_manifest.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(man[0]), delimiter="\t"); w.writeheader(); w.writerows(man)
    # coverage table per type
    cov = {}
    for m in man: cov.setdefault(m["project"], {"T": set(), "N": set()})[m["side"]].add(m["patient"])
    log("  type  pairs_T_RNA  pairs_N_RNA  complete_T+N")
    for t in TYPES13:
        T, N = cov.get(t, {"T": set(), "N": set()}).values()
        log(f"  {t:5s} {len(T):11d} {len(N):12d} {len(T & N):13d}")
    os.makedirs(rna_dir, exist_ok=True)
    failed = []
    for i, m in enumerate(man, 1):
        try:
            st = fetch(f"{API}/data/{m['file_id']}", f"{rna_dir}/{m['file_id']}.tsv", m["md5"], m["size"])
        except SystemExit as e:
            st = "FAILED"; failed.append(m["file_id"]); log(f"  {e}")
        if i % 25 == 0 or i == len(man): log(f"  RNA {i}/{len(man)} ({st})")
    if failed:
        open(f"{outdir}/rna_failed.txt", "w").write("\n".join(failed) + "\n")
        raise SystemExit(f"{len(failed)} RNA files failed (list: {outdir}/rna_failed.txt); re-run to retry, matrix not built")
    # log2(TPM+1) matrix, protein-coding genes, columns = <patient>_<T|N>
    log("  building log2(TPM+1) matrix")
    genes, cols, mat = None, [], []
    for m in man:
        g, v = [], []
        with open(f"{rna_dir}/{m['file_id']}.tsv") as f:
            for line in f:
                if line.startswith(("#", "gene_id", "N_")): continue
                p = line.rstrip("\n").split("\t")  # gene_id gene_name gene_type unstranded ... tpm_unstranded(6)
                if p[2] != "protein_coding": continue
                g.append(f"{p[0]}|{p[1]}"); v.append(p[6])
        if genes is None: genes = g
        elif g != genes: raise SystemExit(f"gene order differs: {m['file_id']}")
        cols.append(f"{m['patient']}_{m['side']}"); mat.append(v)
    import math
    with gzip.open(f"{outdir}/rna_log2tpm_proteincoding.tsv.gz", "wt") as f:
        f.write("gene\t" + "\t".join(cols) + "\n")
        for j, gname in enumerate(genes):
            f.write(gname + "\t" + "\t".join(f"{math.log2(float(c[j]) + 1):.4f}" for c in mat) + "\n")
    log(f"  matrix: {len(genes)} genes x {len(cols)} samples -> {outdir}/rna_log2tpm_proteincoding.tsv.gz")


def block_tables(outdir):
    for name, uuid in GDC_TABLES.items():
        log(f"table {name}: {fetch(f'{API}/data/{uuid}', f'{outdir}/{name}')}")
    for name, url in OTHER_URLS.items():
        try: log(f"table {name}: {fetch(url, f'{outdir}/{name}')}")
        except SystemExit: log(f"  MANUAL: download {url} -> {outdir}/{name}")


def block_genesets(outdir):
    for gs in GENESETS:
        out = f"{outdir}/{gs}.json"
        try: log(f"geneset {gs}: {fetch(MSIGDB.format(gs), out)}")
        except SystemExit: log(f"  MANUAL: {MSIGDB.format(gs)} -> {out}")


def block_clinical(pairs, outdir):
    log("clinical: GDC diagnoses (grade, stage, Gleason) for all patients")
    filt = {"op": "in", "content": {"field": "submitter_id", "value": sorted(pairs)}}
    hits = gdc_query("cases", filt, ["submitter_id", "project.project_id", "diagnoses.tumor_grade",
                     "diagnoses.ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_t",
                     "diagnoses.primary_gleason_grade", "diagnoses.secondary_gleason_grade",
                     "diagnoses.primary_diagnosis", "diagnoses.morphology"])
    keys = ["tumor_grade", "ajcc_pathologic_stage", "ajcc_pathologic_t", "primary_gleason_grade",
            "secondary_gleason_grade", "primary_diagnosis", "morphology"]
    with open(f"{outdir}/gdc_diagnoses.tsv", "w") as f:
        f.write("patient\tproject\t" + "\t".join(keys) + "\n")
        for h in hits:
            d = (h.get("diagnoses") or [{}])[0]
            f.write(f"{h['submitter_id']}\t{h['project']['project_id']}\t" +
                    "\t".join(str(d.get(k, "")) for k in keys) + "\n")
    log(f"  {len(hits)} cases -> {outdir}/gdc_diagnoses.tsv")


def check_coverage(pairs, outdir):
    """Per type: how many of our pairs have ABSOLUTE ploidy/WGD, aneuploidy score, telomere, proliferation score."""
    def ids_in(path, opener=open):
        try:
            with opener(path, "rt") as f: return set(re.findall(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", f.read()))
        except FileNotFoundError: return None
    src = {"ABSOLUTE": ids_in(f"{outdir}/absolute_purity_ploidy.tsv"),
           "aneuploidy": ids_in(f"{outdir}/aneuploidy_arm_calls.tsv"),
           "panimmune": ids_in(f"{outdir}/panimmune_scores160.tsv.gz", gzip.open),
           "subtypes": ids_in(f"{outdir}/tcga_subtypes.tsv")}
    log("coverage of our patients (any sample; tumour-sample match checked at analysis): " +
        ", ".join(f"{k} {len(v & set(pairs)) if v is not None else 'missing'}" for k, v in src.items()))
    sc = f"{outdir}/panimmune_scores160.tsv.gz"
    if os.path.exists(sc):
        with gzip.open(sc, "rt") as f:
            rows = [l.split("\t")[0] for l in f if re.search("prolif", l.split("\t")[0], re.I)]
        log(f"  panimmune proliferation rows: {rows or 'NONE FOUND (use Thorsson 2018 Table S1 by hand)'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs/pairs_gdc.tsv")
    ap.add_argument("--outdir", default="data/c8")
    ap.add_argument("--rna_dir", default="/home/chattopa/data_storage/TCGA_GDC_RNA_STAR")
    ap.add_argument("--only", choices=["rna", "tables", "genesets", "clinical"])
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pairs = read_pairs(a.pairs)
    log(f"{len(pairs)} patients read from {a.pairs} "
        f"({sum(1 for t, n in pairs.values() if t and n)} with both methylation barcodes)")
    if not pairs: sys.exit("no TCGA barcodes found in pairs file")
    if a.only in (None, "tables"): block_tables(a.outdir)
    if a.only in (None, "genesets"): block_genesets(a.outdir)
    if a.only in (None, "clinical"): block_clinical(pairs, a.outdir)
    if a.only in (None, "rna"): block_rna(pairs, a.outdir, a.rna_dir)
    check_coverage(pairs, a.outdir)
    log("done")
