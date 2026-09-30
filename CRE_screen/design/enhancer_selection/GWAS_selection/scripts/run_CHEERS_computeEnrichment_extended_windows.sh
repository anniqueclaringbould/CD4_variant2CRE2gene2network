#Submit python script for CHEERS enrichment
#https://github.com/TrynkaLab/CHEERS

conda activate cheers

#Input info
ROOT=/g/steinmetz/project/otar
DIR=$ROOT/select_enhancer_peaks
SCRIPTS=$ROOT/scripts
GWAS=$ROOT/data/gwas_selection
PEAKS=$ROOT/data/trynka_lab/CHEERS_peak_files/PEAKS

#Make output directory
mkdir $DIR/enrichments

	#Run script for each disease separately
	TRAITS=$(cat $GWAS/disease_names.tsv)

	for TRAIT in $TRAITS

	do

		#Run script for ATAC and H3K27ac peaks separately

			mkdir $DIR/enrichments/ATAC
			mkdir $DIR/enrichments/H3K27ac

			#Run script seperately for extended windows (only for ATAC data)

			for WINDOW in 1000 3000 5000

			do

				mkdir $DIR/enrichments/ATAC/$WINDOW

				echo "
##############################################################
NOW RUNNING: $TRAIT for ATAC peaks (WINDOW of $WINDOW)
##############################################################
				"

					#--input INPUT Text file containing peak coordinates and specificity scores for each of the analyzed samples (default: None)
					#--ld LD Directory with LD information for each SNP (default: None)
					#--snp_list SNP_LIST list of SNPs if CHEERS is used on finemapped set (default: None)
					#--trait TRAIT Name of the analyzed trait (default: None)
					#--outdir OUTDIR Directory where to output results (default: None)

					python $SCRIPTS/CHEERS_computeEnrichment_extended_windows.py \
					--input $DIR/normalised_peaks_across_conditions/ATAC_counts_normToMax_quantileNorm_euclideanNorm_extended_${WINDOW}_filtered.txt \
					--ld $GWAS/$TRAIT/LD/ \
					--trait $TRAIT \
					--outdir $DIR/enrichments/ATAC/$WINDOW/

			done

			echo "
###################################################
NOW RUNNING: $TRAIT for H3K27ac peaks
###################################################
			"

				#--input INPUT Text file containing peak coordinates and specificity scores for each of the analyzed samples (default: None)
				#--ld LD Directory with LD information for each SNP (default: None)
				#--snp_list SNP_LIST list of SNPs if CHEERS is used on finemapped set (default: None)
				#--trait TRAIT Name of the analyzed trait (default: None)
				#--outdir OUTDIR Directory where to output results (default: None)

				python $SCRIPTS/CHEERS_computeEnrichment.py \
				--input $DIR/normalised_peaks_across_conditions/H3K27ac_counts_normToMax_quantileNorm_euclideanNorm_filtered.txt \
				--ld $GWAS/$TRAIT/LD/ \
				--trait $TRAIT \
				--outdir $DIR/enrichments/H3K27ac/

done


#python CHEERS_computeEnrichment.py --input $DIR/normalised_peaks/ATAC_counts_normToMax_quantileNorm_euclideanNorm.txt --ld $GWAS/Asthma/LD --trait Asthma --outdir $DIR/enrichments/
