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
- **04_01_Weinstock_comparison.ipynb**: Compares experimental DEGs to Weinstock et al. reference and SCEPTRE results, generating comparison plots/tables.
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
