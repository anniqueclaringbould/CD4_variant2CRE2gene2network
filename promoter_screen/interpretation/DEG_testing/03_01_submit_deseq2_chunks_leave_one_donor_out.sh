#!/bin/bash
#SBATCH --job-name=deseq2_lodo_chunk_%a
#SBATCH --output=/g/stegle/schrod/code/TCell/Donor/slurm_logs/deseq2_lodo_chunk_%a_%A.out
#SBATCH --error=/g/stegle/schrod/code/TCell/Donor/slurm_logs/deseq2_lodo_chunk_%a_%A.err
#SBATCH --time=4:00:00
#SBATCH --mem=70G
#SBATCH --cpus-per-task=20
#SBATCH --array=0-30   ### Adjust this range based on number of chunks you need


#source the correct environment


# Calculate start and end index for this chunk
N_CHUNKS=5            #Set approptiately!!!!!!!! for 50 ->5, for 250 -> 1
START_P=$(( SLURM_ARRAY_TASK_ID * N_CHUNKS ))
END_P=$(( (SLURM_ARRAY_TASK_ID + 1) * N_CHUNKS ))


# Run the script for this chunk. The left out donor is passed on the sbatch command line via $@
python /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_run_DESeq2_batched.py \
    --chunk-start ${START_P} \
    --chunk-end ${END_P} \
    $@
echo "Completed chunk ${SLURM_ARRAY_TASK_ID}"



# One submission per left out donor, each with its own output directory. Filter settings are the same
# as in the full analysis, since two donors still give up to 6 pseudobulks per target.
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_leave_one_donor_out.sh --exclude-donor A --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Donor/pseudobulk_deseq2_without_donor_A
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_leave_one_donor_out.sh --exclude-donor B --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Donor/pseudobulk_deseq2_without_donor_B
#sbatch /g/stegle/schrod/code/CD4_variant2CRE2gene2network/promoter_screen/interpretation/DEG_testing/03_01_submit_deseq2_chunks_leave_one_donor_out.sh --exclude-donor C --chunk-size 50 --min-bulks 3 --ribosomal --output-dir /g/stegle/schrod/code/TCell/Donor/pseudobulk_deseq2_without_donor_C
