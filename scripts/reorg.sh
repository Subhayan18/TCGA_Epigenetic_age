#!/usr/bin/env bash
# Reorganise Meth_TCGA_GSE213478/ into data/ scripts/ results/. Dry run by default; "bash reorg.sh go" to move.
# Run from Meth_TCGA_GSE213478/ AFTER the c1_phylo_units.py clusters run has finished.
# Leaves top-level symlinks for paths that older scripts hard-code.
set -u; GO=${1:-}
mv_to() { local d=$1; shift; for f in "$@"; do [ -e "$f" ] && [ ! -L "$f" ] || continue
  echo "$f -> $d/"; [ "$GO" = go ] && mkdir -p "$d" && mv -n "$f" "$d/"; done; }
mv_to scripts                *.py *.sh
mv_to data/pairs             pairs_rebuilt.tsv pairs_gdc.tsv gdc_aliquot_map.tsv gdc_barcodes.tsv excluded_legacy_pairs.tsv
mv_to data/reference         c2_ref
mv_to results/c0_cohort      table1_rebuilt.tsv Supp_cohort_flow.png Supp_excluded_pairs.csv Supp_types_below_threshold.csv
mv_to results/c1_clocks      paired_results_rebuilt.csv pair_table_corrected.csv accel_stats_corrected.csv Figure[123]_corrected.png \
                             c1_out c1_out_v2 Table_mitotic_clocks_tumour_vs_normal.csv Table_normal_mitotic_clock_vs_age.csv Fig_Horvath_vs_mitotic_delta.*
mv_to results/c2_composition c2_out Table_clock_composition_adjusted.csv Supp_composition_by_type.csv check_Epi_vs_CPE.csv
mv_to results/c4_xchrom      phylo_out/C4_2_*
mv_to results/c3_phylo       phylo_out
mv_to archive                c2_bundle.tar.gz
mv_to legacy_student         methylation_betas Adult_DNAmAge_figures_data.xlsx DNAmAge_with_clinical.csv
[ "$GO" = go ] || { echo "(dry run: nothing moved; rerun with: bash reorg.sh go)"; exit 0; }
[ -d results/c3_phylo/phylo_out ] && mv -n results/c3_phylo/phylo_out/* results/c3_phylo/ && rmdir results/c3_phylo/phylo_out
for l in data/pairs/pairs_gdc.tsv results/c1_clocks/pair_table_corrected.csv results/c1_clocks/c1_out \
         results/c2_composition/c2_out data/reference/c2_ref; do [ -e "$l" ] && ln -sfn "$l" "$(basename "$l")"; done
[ -d results/c3_phylo ] && ln -sfn results/c3_phylo phylo_out
echo "done; leftovers at top level:"; ls -p | grep -v /
