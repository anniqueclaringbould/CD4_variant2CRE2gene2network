#!/bin/bash
#SBATCH --job-name=batchjob_ss
#SBATCH --output=slurm_logs/batchjob_ss_%A_%a.out
#SBATCH --error=slurm_logs/batchjob_ss_%A_%a.err
#SBATCH --array=0-9                 # We have 7669 perturbations # 0-9 for 800 batches
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem


# Calculate start and end index for this job
BATCH_SIZE=800   #TODO: change back to 800   # 17 minutes for 50 targets!
START_P=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
END_P=$(( (SLURM_ARRAY_TASK_ID + 1) * BATCH_SIZE ))

python 01_02_run_Wilcoxon.py --start-idx $START_P --end-idx $END_P $@





# All cells:
# sbatch 01_02_run_single_cell_Wilcoxon.sh --adata-path /g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg.h5ad --results-dir /g/stegle/schrod/code/TCell/results_single_cell/Wilcoxon_single_cell_final --cutoff 5
