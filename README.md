# TCGA_Epigenetic_age

Code, small input tables and result tables for **Mukherji & Chattopadhyay** (manuscript in revision): epigenetic and mitotic ageing in paired tumour and adjacent normal tissue across 20 TCGA cancer types.

## Overview

The study asks three questions:

1. **Do tumours look epigenetically older than the normal tissue next to them?** Horvath DNAmAge and four mitotic clocks are compared between each tumour and its own adjacent normal sample.
2. **Is the tumour–normal difference real, or produced by composition, purity, heterogeneity or proliferation?** Cell-type composition (EpiDISH), tumour purity (CPE), per-patient heterogeneity of methylation change, RNA proliferation scores, ploidy/WGD, aneuploidy, telomere length, grade, stage and subtype are each tested as explanations.
3. **Is methylation change on the X chromosome enriched, and does that depend on sex or X copy number?**

A paediatric extension (C1.5) applies the same mitotic clocks to paediatric tumours, including relapse pairs.

### Cohort

| | |
|---|---|
| Platform | Illumina HumanMethylation450 (TCGA, GDC open access) |
| Design | Paired primary tumour (`-01`) vs solid tissue normal (`-11`), one pair per patient |
| Pairs | 531 across 20 cancer types (`results/c0_cohort/table1_rebuilt.tsv`) |
| Main analysis | 506 pairs in the 13 types with ≥10 pairs |
| Excluded | 14 legacy pairs that could not be verified against GDC (`data/pairs/excluded_legacy_pairs.tsv`) |

## Repository layout

```
.
├── c1/                     Mitotic clock implementation + clock CpG lists
│   ├── mitotic_clocks.py
│   └── clock_cpgs/         epiTOC1, epiTOC2 (coefficients), HypoClock, stemTOC
├── Horvath1.csv            Horvath (2013) 353-CpG coefficients
├── scripts/                Analysis pipeline, prefixed by manuscript section (c1_ … c15_)
├── data/
│   ├── pairs/              Tumour–normal pair tables, GDC aliquot/file maps
│   ├── reference/c2_ref/   EpiDISH reference centroids
│   ├── gdc_clinical/       Case-level clinical data (GDC)
│   ├── gdc_ascat/          ASCAT segment manifest
│   └── c8/                 PanCanAtlas/CDR tables, subtypes, Hallmark gene sets, RNA manifest
├── results/
│   ├── c0_cohort/          Cohort table
│   ├── c1_clocks/          Clock estimates, tumour vs normal tests, paediatric clocks
│   ├── c2_composition/     Composition, purity, composition-adjusted clocks
│   ├── c3_phylo/           Heterogeneity of methylation change (500 bp clusters, 5 Mb windows)
│   ├── c4_xchrom/          X-chromosome enrichment
│   ├── c8_prolif/          Proliferation, genomic, telomere, clinical models
│   └── figures/            Main figures
├── archive/, legacy_student/, Phylogenetics.r
│                           Earlier R cladogram analyses (superseded; kept for provenance)
└── requirements.txt
```

Every `Supp_Tables_C*.xlsx` workbook has one sheet per table, each with a bold title row and a note row.

## Pipeline

All scripts are run from the repository root (originally `Meth_TCGA_GSE213478/`). Run `python3 <script> -h` to see a script's arguments. Scripts named `*_supp_tables.py` and `c8_10_outputs.py` only read existing results and never refit models.

### Cohort

`data/pairs/pairs_rebuilt.tsv` is the verified pair table. It was rebuilt from GDC so that every tumour and normal aliquot is checked against its barcode and file UUID (`gdc_barcodes.tsv`, `gdc_aliquot_map.tsv`).

### Clocks

| Script | What it does | Main output |
|---|---|---|
| `scripts/horvath.py` | Horvath DNAmAge per sample: weighted sum + intercept, then anti-log-linear transform; a missing CpG is imputed with the across-sample mean | `paired_results_rebuilt.csv` |
| `c1/mitotic_clocks.py` | Python port of [EpiMitClocks](https://github.com/aet21/EpiMitClocks): **epiTOC1** (mean β, 385 PCGT CpGs), **epiTOC2** (tnsc, 163 CpGs; irS = tnsc/age), **HypoClock** (1 − mean β, 678 solo-WCGW CpGs), **stemTOC** (95th percentile β, 371 CpGs) | `c1_out/mitotic_clocks_{samples,pairs}.tsv`, `clock_coverage.tsv` |
| `scripts/c1_tables.py` | Tumour vs normal tests per type; normal-tissue clock vs chronological age | `Table_mitotic_clocks_tumour_vs_normal.csv`, `Table_normal_mitotic_clock_vs_age.csv`, `accel_stats_corrected.csv` |
| `scripts/c1_supp_tables.py` | Assembles the C1 supplementary workbook | `Supp_Tables_C1.xlsx` |

### Heterogeneity of methylation change

Each patient's probes are grouped into units: probe clusters with gaps ≤500 bp, or fixed 5 Mb windows. For each unit, Δ = mean β(T) − mean β(N).

| Script | What it does |
|---|---|
| `scripts/c1_hetero.py` | Per-patient metrics at \|Δ\| thresholds of 0.05/0.10/0.20: fraction changed, discordance, relative change (Δ scaled to the largest change possible), q95 of relative change, subclonal fraction, k by ICL on a 1-D Gaussian mixture. Output: `hetero_<mode>_patients.tsv` |
| `scripts/c1_hetero_v1.py` | Earlier version of `c1_hetero.py`; its output feeds `c3_supp_tables.py --h500v1` |
| `scripts/c1_phylo_units.py` | Independent re-implementation of the thesis "cladogram" (trunk/private) method. `--mode windows`: UPGMA, k by silhouette. `--mode clusters`: GMM, k by BIC |
| `scripts/c3_hetero_models.py` | metric ~ Δclock + mean/rms/sd Δβ + CPE + age + sex [+ type FE]; HC3 standard errors, BH correction (Supp. Table 23) |
| `scripts/c3_5_orientation.py` | Per-type summary of the 5 Mb replica, compared with the original thesis numbers |
| `scripts/c3_supp_tables.py` | `Supp_Tables_C3.xlsx` |

### Composition and purity

| Script | What it does |
|---|---|
| `scripts/c2_composition.py` | EpiDISH RPC (Python port: Huber robust regression on reference centroids). Uses Epi/Fib/IC for all types and Epi/Fib/Fat/IC for BRCA; CPE purity from Aran 2015. Composition-adjusted clock: Δclock ~ ΔFib + ΔIC [+ ΔFat], where the intercept is the T−N difference at equal composition |
| `scripts/c2_4_hetero_purity.py` | Tests whether the heterogeneity metrics are driven by purity or composition (per type: Spearman; pooled: model with type FE). LIHC is excluded from the composition terms |
| `scripts/c2_supp_tables.py` | `Supp_Tables_C2.xlsx` |

### X chromosome

| Script | What it does |
|---|---|
| `scripts/c4_fetch_gdc.py` | Downloads sex and ASCAT allele-specific copy number from GDC; summarises X copy-number state per patient |
| `scripts/c4_2_x_enrichment.py` | X event rate divided by autosomal rate (events per tested unit); paired Wilcoxon within each type, BH across types |
| `scripts/c4_1_3_x_sex_cna.py` | X enrichment stratified by sex and by X copy-number state (strata with n ≥ 10) |
| `scripts/c4_supp_tables.py` | `Supp_Tables_C4.xlsx` |

### What does the mitotic-clock advance track?

| Script | What it does |
|---|---|
| `scripts/c8_1_download.py` | Resumable download of RNA-seq, PanCanAtlas tables, gene sets and clinical data (`--only rna\|tables\|genesets\|clinical`) |
| `scripts/c8_2_prolif_scores.py` | Proliferation scores: core (MKI67, PCNA, MCM2–7), Yang 2016 9-gene mitotic index, Hallmark E2F/G2M (singscore-style) |
| `scripts/c8_3_genomic.py` | ABSOLUTE ploidy/WGD, aneuploidy score, PanCanAtlas QC flags |
| `scripts/c8_4_telomere.py` | Telomere length, log2(TL tumour / TL normal), from Barthel 2017 |
| `scripts/c8_5_grade_subtype.py` | Availability of grade, stage (TCGA-CDR, GDC) and mRNA subtype |
| `scripts/c8_6_models.py` | Δclock ~ proliferation + purity + composition + Δβ + age + sex; incremental value for grade, stage, WGD and aneuploidy (LR tests); block contributions. `--check` only checks which input columns matched |
| `scripts/c8_10_outputs.py` | Figure 6 composite (SVG/PDF) and `Supp_Tables_C8.xlsx` |

### Paediatric tumours

| Script | What it does |
|---|---|
| `scripts/c15_ped_cpgmatrix.py` | β values at the union of all clock CpGs per sample, so clocks can be recomputed on a CpG set shared by 450K, EPIC and EPICv2 |
| `scripts/c15_ped_clocks.py` | Mitotic clocks per paediatric sample; results are in `results/c1_clocks/pediatric/` |

### Utility

`scripts/reorg.sh` records how the original working directory was restructured into `data/ scripts/ results/`. It is a dry run unless called with `go`.

## Data not included

Raw data are not distributed. To reproduce from scratch, download:

- TCGA 450K β-value files, clinical data and ASCAT segments from the [GDC portal](https://portal.gdc.cancer.gov/) (open access). File UUIDs are listed in `data/pairs/gdc_barcodes.tsv`.
- RNA-seq and PanCanAtlas tables with `scripts/c8_1_download.py`.
- Paediatric methylation data from GEO [GSE213478](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213478).

## Requirements

- Python ≥ 3.9: `pip install -r requirements.txt` (numpy, pandas, scipy, statsmodels, matplotlib, openpyxl)
- Legacy R scripts only: tidyverse, ape, phangorn, ggdendro, cowplot, openxlsx

## Methods references

- Horvath S. *Genome Biol* 2013;14:R115 (DNAmAge)
- Yang Z et al. *Genome Biol* 2016;17:205 (epiTOC1, mitotic index)
- Teschendorff AE. *Genome Med* 2020;12:56 (epiTOC2, HypoClock)
- Zhu T et al. *Nat Commun* 2024;15:4211 (stemTOC)
- Teschendorff AE et al. *BMC Bioinformatics* 2017;18:105 (EpiDISH)
- Aran D et al. *Nat Commun* 2015;6:8971 (CPE purity)
- Taylor AM et al. *Cancer Cell* 2018;33:676 (aneuploidy)
- Barthel FP et al. *Nat Genet* 2017;49:349 (telomere length)
- Liu J et al. *Cell* 2018;173:400 (TCGA-CDR)
- Foroutan M et al. *BMC Bioinformatics* 2018;19:404 (singscore)

## Citation

Mukherji A, Chattopadhyay S. Manuscript in revision. A citation will be added on publication.
