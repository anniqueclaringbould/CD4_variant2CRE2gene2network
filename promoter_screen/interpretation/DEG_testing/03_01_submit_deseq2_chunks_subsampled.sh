#!/bin/bash
#SBATCH --job-name=deseq2_sub_chunk_%a
#SBATCH --output=slurm_logs/deseq2_sub_chunk_%a_%A.out
#SBATCH --error=slurm_logs/deseq2_sub_chunk_%a_%A.err
#SBATCH --time=4:00:00
#SBATCH --mem=70G
#SBATCH --cpus-per-task=20
#SBATCH --array=0-30   ### Adjust this range based on number of chunks you need


#source the correct environment


# Calculate start and end index for this chunk
N_CHUNKS=5            #Set approptiately!!!!!!!! for 50 ->5, for 250 -> 1
START_P=$(( SLURM_ARRAY_TASK_ID * N_CHUNKS ))
END_P=$(( (SLURM_ARRAY_TASK_ID + 1) * N_CHUNKS ))


# Run the script for this chunk. The kept cell fraction is passed on the sbatch command line via $@
python /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_run_DESeq2_batched.py \
    --chunk-start ${START_P} \
    --chunk-end ${END_P} \
    $@
echo "Completed chunk ${SLURM_ARRAY_TASK_ID}"



# Random removal of 10%, 25% and 50% of the cells, all three donors kept. Same filter settings as the
# full analysis, so the results are directly comparable to pseudobulk_deseq2_chunks50_min3bulks_ribosomal.
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.90 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep90_seed0
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.75 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep75_seed0
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.50 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep50_seed0

# Fractions matched to what each leave-one-donor-out run actually retained (79.4%, 65.0%, 55.6% of
# cells). These are the cells-matched controls that separate loss of power from donor identity.
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.794 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep79_seed0
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.650 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep65_seed0
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.556 --seed 0 --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep56_seed0

# Repeat a fraction with different seeds to get a spread rather than a single draw:
#for S in 0 1 2; do sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_subsampled.sh --keep-frac 0.50 --seed $S --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Subsample/pseudobulk_deseq2_keep50_seed$S; done

# Remember to create the slurm_logs directory in the submission directory first, SLURM does not.
