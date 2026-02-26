#!/bin/bash
#SBATCH --job-name=NT_split_%a
#SBATCH --output=slurm_logs/NT_split_%a_%A.out
#SBATCH --error=slurm_logs/NT_split_%a_%A.err
#SBATCH --time=00:15:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=20
#SBATCH --array=0-299     #00  # Run 300 iterations (adjust as needed)


# Source your environment (adjust as needed)
# source activate your_environment

# Run the script for this iteration
python /g/stegle/schrod/code/TCell/03_02_DESeq2_calibration.py \
    --iteration ${SLURM_ARRAY_TASK_ID} \
    --n-guides 3 \
    --n-cpus-fit 20 \
    --n-cpus-test 20 \
    --output-dir /g/stegle/schrod/code/TCell/pseudobulk_deseq2_NT_ribosomal

echo "Completed iteration ${SLURM_ARRAY_TASK_ID}"


# sbatch 03_02_run_DESeq2_calibration.sh