# Additional_DEG_methods README

This folder contains scripts, notebooks, and batch scripts for running and evaluating alternative differential expression (DEG) methods, primarily based on Wilcoxon tests, for single-cell and pseudobulk data.

## Required Input Files

- **promoter_single_guide_hvg_raw_counts.h5ad**: Preprocessed AnnData for pseudobulk construction
- **promoter_single_guide_hvg.h5ad**: Preprocessed AnnData for downstream analyses


## File Descriptions

- **01_01_Wilcoxon_single_cell_calibration.py**: Performs calibration of single-cell Wilcoxon DEG results.
- **01_02_run_single_cell_Wilcoxon.sh**: Batch script to run single-cell Wilcoxon DEG analysis.
- **01_02_run_Wilcoxon.py**: Runs Wilcoxon DEG analysis on single-cell data.
- **01_03_evaluate_calibration_Wilcoxon.ipynb**: Evaluates and visualizes calibration of single-cell Wilcoxon results.
- **02_01_build_pseudobulks_for_Wilcoxon.ipynb**: Prepares pseudobulk data for Wilcoxon DEG analysis.
- **02_02_run_pseudobulk_Wilcoxon.sh**: Batch script to run Wilcoxon DEG analysis on pseudobulk data.
- **02_03_run_Wilcoxon_bulk_calibration_test.sh**: Batch script to run calibration tests for Wilcoxon on bulk/pseudobulk data.
- **02_03_Wilcoxon_bulk_calibration.py**: Performs calibration of Wilcoxon DEG results on pseudobulk data.
- **02_04_evaluate_pseudobulk_calibration.ipynb**: Evaluates and visualizes calibration of pseudobulk Wilcoxon results.
- **03_00_compare_different_DEG_methods.ipynb**: Compares results from different DEG methods (single-cell and pseudobulk).

---
