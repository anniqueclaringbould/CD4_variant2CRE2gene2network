#!/bin/bash
#SCRIPT TO MAKE AND SUBMIT SCRIPTS TO PROCESS SEQUENCING DATA BY BD RHAPSODY PIPELINE

#All samples that have A and B side (no P7)

#add sample names for files without P addition
samples=(A10_N B10_N C10_N)

# Associate expected cell numbers with each sample
declare -A cell_counts=( ["A1"]=123647 ["A2"]=39662 ["A3"]=54709 ["A4"]=53343 ["A5"]=40250 ["A6"]=52562 ["A7"]=54513 ["A8_N"]=52465 ["A9_N"]=46081 ["A10_N"]=41962 ["A11_N"]=112908 ["A12"]=93796 ["A13"]=100884 ["A14"]=93705 ["A15"]=93618 ["A16"]=88975 ["B1"]=51666 ["B2"]=48399 ["B3"]=48832 ["B4"]=47326 ["B5"]=52612 ["B6"]=51224 ["B7"]=51445 ["B8_N"]=113410 ["B9_N"]=43136 ["B10_N"]=50029 ["B11_N"]=41468 ["B12"]=47659 ["B13"]=51619 ["B14"]=48863 ["B15"]=48658 ["B16"]=46851 ["C1"]=49901 ["C2"]=50124 ["C3"]=45854 ["C4"]=51978 ["C5"]=56722 ["C6"]=49749 ["C7"]=43981 ["C8_N"]=52684 ["C9_N"]=52915 ["C10_N"]=48082 ["C11_N"]=53320 ["C12"]=49357 ["C13"]=52043 ["C14"]=49901 ["C15"]=48990 ["C16"]=48788 )

#save root folder
ROOT=/g/steinmetz/project/otar


for sample in "${samples[@]}"

do
	echo "Currently making and submitting scripts for ${sample}"
	echo "Expected cell count: ${cell_counts[$sample]}"

	#make scripts per sample
	echo "#!/bin/bash
#SBATCH -A lsteinme            		 	# group to which you belong
#SBATCH -J ${sample}					# job_name
#SBATCH -n 32                		  	# number of cores
#SBATCH -N 1                 		    # number of nodes
#SBATCH -p bigmem 
#SBATCH -C naples
#SBATCH --mem 128312            		# memory pool for all cores
#SBATCH -t 3-00:00:00           		# runtime limit (D-HH:MM:SS)
#SBATCH -o $ROOT/software/BD2.2/full_data_per_panel/logs/${sample}.out
#SBATCH -e $ROOT/software/BD2.2/full_data_per_panel/logs/${sample}.err
#SBATCH --mail-type=END,FAIL        	# notifications for job done & fail
#SBATCH --mail-user=dewi.moonen@embl.de # send-to address

$ROOT/software/rhapsodyPipeline-2.2/rhapsody pipeline \
--outdir $ROOT/software/BD2.2/full_data_per_panel/output/${sample}/ \
$ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_${sample}.yml" > $ROOT/software/BD2.2/full_data_per_panel/scripts/rhapsody_pipeline2.2_Enhancerscreen_${sample}.sh
	
	#make yml inputs per sample
    sed -e "s/SAMPLE/${sample}/g" -e "s/EXPECTED_CELLS/${cell_counts[$sample]}/g" $ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_template.yml > $ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_${sample}.yml

	#submit script per sample
	sbatch $ROOT/software/BD2.2/full_data_per_panel/scripts/rhapsody_pipeline2.2_Enhancerscreen_${sample}.sh

done

#All samples that have A and B side + P7

#add sample names for files with P addition
samples=(A10_N)

#save root folder
ROOT=/g/steinmetz/project/otar


for sample in "${samples[@]}"

do
	echo "Currently making and submitting scripts for ${sample}"
	echo "Expected cell count: ${cell_counts[$sample]}"

	#make scripts per sample
	echo "#!/bin/bash
#SBATCH -A lsteinme            		 	# group to which you belong
#SBATCH -J ${sample}					# job_name
#SBATCH -n 32                		  	# number of cores
#SBATCH -N 1                 		    # number of nodes
#SBATCH -p bigmem 
#SBATCH -C naples
#SBATCH --mem 128312            		# memory pool for all cores
#SBATCH -t 3-00:00:00           		# runtime limit (D-HH:MM:SS)
#SBATCH -o $ROOT/software/BD2.2/full_data_per_panel/logs/${sample}.out
#SBATCH -e $ROOT/software/BD2.2/full_data_per_panel/logs/${sample}.err
#SBATCH --mail-type=END,FAIL        	# notifications for job done & fail
#SBATCH --mail-user=dewi.moonen@embl.de # send-to address

$ROOT/software/rhapsodyPipeline-2.2/rhapsody pipeline \
--outdir $ROOT/software/BD2.2/output/${sample}/ \
$ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_${sample}.yml" > $ROOT/software/BD2.2/full_data_per_panel/scripts/rhapsody_pipeline2.2_Enhancerscreen_${sample}.sh
	
	#make yml inputs per sample
    sed -e "s/SAMPLE/${sample}/g" -e "s/EXPECTED_CELLS/${cell_counts[$sample]}/g" $ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_template_withP.yml > $ROOT/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_${sample}.yml

	#submit script per sample
	sbatch $ROOT/software/BD2.2/full_data_per_panel/scripts/rhapsody_pipeline2.2_Enhancerscreen_${sample}.sh

done