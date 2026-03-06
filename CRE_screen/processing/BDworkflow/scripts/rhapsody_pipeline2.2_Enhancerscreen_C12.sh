#!/bin/bash
#SBATCH -A lsteinme            		 	# group to which you belong
#SBATCH -J C12					# job_name
#SBATCH -n 32                		  	# number of cores
#SBATCH -N 1                 		    # number of nodes
#SBATCH -p bigmem 
#SBATCH -C naples
#SBATCH --mem 128312            		# memory pool for all cores
#SBATCH -t 3-00:00:00           		# runtime limit (D-HH:MM:SS)
#SBATCH -o /g/steinmetz/project/otar/software/BD2.2/full_data_per_panel/logs/C12.out
#SBATCH -e /g/steinmetz/project/otar/software/BD2.2/full_data_per_panel/logs/C12.err
#SBATCH --mail-type=END,FAIL        	# notifications for job done & fail
#SBATCH --mail-user=dewi.moonen@embl.de # send-to address

/g/steinmetz/project/otar/software/rhapsodyPipeline-2.2/rhapsody pipeline --outdir /g/steinmetz/project/otar/software/BD2.2/full_data_per_panel/output/C12/ /g/steinmetz/project/otar/software/BD2.2/full_data_per_panel/input/pipeline_inputs_2.2_Enhancerscreen_C12.yml
