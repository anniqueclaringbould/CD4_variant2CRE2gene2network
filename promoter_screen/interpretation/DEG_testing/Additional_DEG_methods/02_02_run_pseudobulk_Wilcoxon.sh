#!/bin/bash
#SBATCH --job-name=batchjob_pb
#SBATCH --output=slurm_logs/batchjob_pseudobulk_%A_%a.out
#SBATCH --error=slurm_logs/batchjob_pseudobulk_%A_%a.err
#SBATCH --array=0-9                 # We have 7669 perturbations # 0-9 for 800 batches
#SBATCH --mem=150G
#SBATCH --time=4:00:00 #24:00:00
#SBATCH --cpus-per-task=1


# Calculate start and end index for this job
BATCH_SIZE=800   #TODO: change back to 800
START_P=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
END_P=$(( (SLURM_ARRAY_TASK_ID + 1) * BATCH_SIZE ))

python 01_02_run_Wilcoxon.py --start-idx $START_P --end-idx $END_P $@


# Rerun pseudobulk:
# sbatch 02_02_run_pseudobulk_Wilcoxon.sh --adata-path /g/stegle/schrod/data/T_Cell/promoter_bulk_profiles_1guide.h5ad --results-dir /g/stegle/schrod/code/TCell/results_pseudobulk/Wilcoxon_bulk_final_cutoff_2 --cutoff 2


