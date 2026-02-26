#!/bin/bash
#SBATCH --job-name=deseq2_chunk_%a
#SBATCH --output=slurm_logs/deseq2_chunk_%a_%A.out
#SBATCH --error=slurm_logs/deseq2_chunk_%a_%A.err
#SBATCH --time=4:00:00
#SBATCH --array=0-19  # Adjust this range based on number of chunks you need
#SBATCH --mem=70G
#SBATCH --cpus-per-task=20
#SBATCH --array=0-30  ### Adjust this range based on number of chunks you need


#source the correct environment


# Calculate start and end index for this chunk
N_CHUNKS=5            #Set approptiately!!!!!!!! for 50 ->5, for 250 -> 1
START_P=$(( SLURM_ARRAY_TASK_ID * N_CHUNKS ))
END_P=$(( (SLURM_ARRAY_TASK_ID + 1) * N_CHUNKS ))


# Run the script for this chunk
python /g/stegle/schrod/code/TCell/07_02_pseudobulk_DESeq2_chunked_parallel.py \
    --chunk-start ${START_P} \
    --chunk-end ${END_P} \
    $@
echo "Completed chunk ${SLURM_ARRAY_TASK_ID}"



# Without ribosomal percentage correction:
#sbatch 03_01_submit_deseq2_chunks.sh --chunk-size 50 --output-dir /g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50_min3bulks --min-bulks 3

# With ribosomal percentage correction:
# sbatch 03_01_submit_deseq2_chunks.sh --chunk-size 50 --output-dir /g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50_min3bulks_ribosomal --min-bulks 3 --ribosomal
