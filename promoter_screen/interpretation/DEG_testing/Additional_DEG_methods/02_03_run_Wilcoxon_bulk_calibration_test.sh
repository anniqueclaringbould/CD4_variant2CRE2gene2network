#!/bin/bash
#SBATCH --job-name=NT_wilcox_%a
#SBATCH --output=slurm_logs/NT_wilcox_%a_%A.out
#SBATCH --error=slurm_logs/NT_wilcox_%a_%A.err
#SBATCH --time=00:15:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-300

# Run the script for this iteration
python 02_03_Wilcoxon_bulk_calibration.py \
    --iteration ${SLURM_ARRAY_TASK_ID} \
    --n-guides 3 \
    --output-dir /g/stegle/schrod/code/TCell/wilcoxon_batch_calibration_NT

echo "Completed iteration ${SLURM_ARRAY_TASK_ID}"

# Calibration for batches
# sbatch 02_03_run_Wilcoxon_bulk_calibration_test.sh
