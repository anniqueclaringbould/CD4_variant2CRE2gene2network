# DEG_testing README

This folder contains scripts and notebooks for preprocessing, pseudobulk construction, differential expression (DE) analysis, DE calibration, and evaluation of knockdown efficiency in the promoter screen.

## Required Datasets and Files
- **combined_anndata.h5ad**: Raw combined single-cell data
- **guide_assignments/matrix/grna_assignment_matrix_colnames.txt**: Guide assignment matrices
- **guide_assignments/matrix/grna_assignment_matrix_rownames.txt**: Guide assignment matrices
- **Weinstock_S4.xlsx**: Reference DEGs from Weinstock et al.

## File Descriptions

- **01_01_preprocess.ipynb**: Preprocesses raw single-cell data and guide assignments, producing pseudobulk-ready AnnData files.
- **01_02_build_pseudobulks.py**: Constructs pseudobulk AnnData objects needed for DESeq2 and generates summary statistics and figures.
- **02_01_Evaluate_knockdown_efficiancy_per_guide.ipynb**: Evaluates knockdown efficiency for each guide and outputs summary tables/plots.
- **03_01_run_DESeq2_batched.py**: Runs DESeq2 in batch mode on pseudobulk data and saves differential expression results.
- **03_02_DESeq2_calibration.py**: Performs calibration of DESeq2 results and saves calibration outputs.
- **03_03_evaluate_DESeq2_calibration.ipynb**: Analyzes and visualizes calibration results from DESeq2.
- **03_04_Evaluate_knockdown_efficiency_per_target.ipynb**: Summarizes knockdown efficiency per target gene using AnnData or DESeq2 results.
- **01_03_add_cell_scores_to_pseudobulks.py**: Aggregates the per-cell cell-state AUC scores to the pseudobulk level so they can be used as covariates in DESeq2.
- **01_04_cell_score_distributions.py**: Distributions of the cell-state AUC scores per cell and per pseudobulk, and how much of each score is explained by donor.
- **03_01_submit_deseq2_chunks_leave_one_donor_out.sh**: Slurm wrapper running the DESeq2 chunks with one donor excluded.
- **03_01_submit_deseq2_chunks_subsampled.sh**: Slurm wrapper running the DESeq2 chunks on random cell subsamples, all donors retained.
- **03_05_compare_leave_one_donor_out.py**: Compares the edges detected in the leave-one-donor-out runs; writes the caches the later donor scripts reuse.
- **03_06_compare_lodo_to_full_run.py**: Compares each leave-one-donor-out run against the full three-donor run and writes the log-fold-change matrix.
- **03_07_donor_contribution.py**: Cells and pseudobulks contributed by each donor.
- **03_08_subsample_summary.py**: Summarises the random cell-subsampling runs on the same universe as the leave-one-donor-out numbers.
- **03_09_compare_subsample_to_full_run.py**: As 03_06, for the random subsampling runs.
- **03_10_cells_vs_discoveries.py**: Cells contributed by each donor against links discovered when that donor is left out, with the matched random baseline.
- **03_11_tyk2_dexi_lfc_across_runs.py**: Log fold changes and standard errors of the TYK2 and DEXI cascade genes across every run.
- **03_12_effect_vs_rank_stability.py**: Agreement of the published edges' effect estimates with the full run, as Pearson r, sign concordance and AUCC.
- **03_13_node_effect_significance_heatmap.py**: Per-gene effects and significance across all runs for the DEXI, SSB and TYK2 downstream sets.
- **03_14_twohop_expression_matched_null.py**: Expression-matched random-seed null for the two-hop network and its programme enrichment.
- **03_15_twohop_null_halfviolin.py**: The two-hop null result as bars against the null's density.
- **04_01_Weinstock_comparison.ipynb**: Compares experimental DEGs to Weinstock et al. reference and SCEPTRE results, generating comparison plots/tables.
- **04_02_Figure4E.py**: Compact TF enrichment figure with motif logos (Figure 4E).
- **utils.py**: Utility functions for normalization, preprocessing, plotting, and DEG result loading.

## Additional differential gene expression methods

Code for additional differential gene expression testing methods explored in the manuscript can be 
found in `Additional_DEG_methods`:

#### SCEPTRE
Scripts to run model calibration checks and the full differential expression analysis using the 
SCEPTRE Nextflow pipeline. Because of the high number of tests, running the pipeline using
large-scale parallization, for instance on a HPC cluster with Slurm as executor, is required. See
the [SCEPTRE Nextflow pipeline](https://timothy-barry.github.io/sceptre-book/at-scale.html) and [Nextflow](https://www.nextflow.io/docs/latest/) documentations for more information.

---
