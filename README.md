# TCGA_Epigenetic_age

Epigenetic (Horvath) and mitotic clocks (epiTOC1/2, HypoClock, stemTOC) in paired tumour–normal TCGA 450K methylation across 20 cancer types, with composition, heterogeneity, X-chromosome and proliferation analyses. Companion code to Mukherji & Chattopadhyay (in revision).

## Layout
| Path | Content |
|---|---|
| `scripts/` | Analysis pipeline, numbered by manuscript section (`c1_`…`c15_`) |
| `c1/` | Mitotic clock implementation (Python port of [EpiMitClocks](https://github.com/aet21/EpiMitClocks)) + clock CpG lists |
| `Horvath1.csv` | Horvath (2013) clock coefficients |
| `data/pairs/` | Tumour–normal pair tables, GDC aliquot/barcode maps |
| `data/reference/` | EpiDISH reference centroids |
| `data/c8/`, `data/gdc_*` | Small GDC/PanCanAtlas-derived inputs (gene sets, clinical, ASCAT) |
| `results/c0`–`c8` | Result tables and supplementary tables per section |
| `results/figures/` | Manuscript figures |
| `archive/`, `legacy_student/`, `Phylogenetics.r` | Earlier R cladogram analyses (not used in the revised manuscript) |

## Pipeline (run from repo root)
| Step | Script | Purpose |
|---|---|---|
| C1.1 | `c1/mitotic_clocks.py`, `scripts/horvath.py` | Clock estimates per sample |
| C1 | `c1_tables.py`, `c1_supp_tables.py` | Tumour vs normal clock tables |
| C1.4 | `c1_hetero.py`, `c1_phylo_units.py` | Per-patient heterogeneity of methylation change |
| C2 | `c2_composition.py`, `c2_4_hetero_purity.py`, `c2_supp_tables.py` | EpiDISH RPC composition, purity-adjusted clocks |
| C3 | `c3_hetero_models.py`, `c3_5_orientation.py`, `c3_supp_tables.py` | Clock ~ heterogeneity models (500 bp vs 5 Mb) |
| C4 | `c4_fetch_gdc.py`, `c4_2_x_enrichment.py`, `c4_1_3_x_sex_cna.py`, `c4_supp_tables.py` | X-chromosome event enrichment by sex/CN |
| C8 | `c8_1_download.py` … `c8_10_outputs.py` | Proliferation, genomic, telomere, grade/subtype models |
| C15 | `c15_ped_cpgmatrix.py`, `c15_ped_clocks.py` | Paediatric tumour clocks |

`scripts/reorg.sh` documents how the original working directory was restructured.

## Data
Raw data are not included. Sources: TCGA via [GDC](https://portal.gdc.cancer.gov/) (450K beta values, clinical, ASCAT), PanCanAtlas supplementary tables, and GEO [GSE213478](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213478). `c4_fetch_gdc.py` and `c8_1_download.py` fetch the inputs they need.

## Requirements
Python ≥3.9: `pip install -r requirements.txt`. Legacy R scripts: tidyverse, ape, phangorn, ggdendro, cowplot, openxlsx.

## Note
Some scripts still contain absolute paths from the original server (e.g. `scripts/horvath.py`); edit before running.
