#!/bin/bash

## Run SCEPTRE Nextflow pipeline to assign gRNAs using the mixture method. Requires that Nextflow
## and R are available in the PATH, and the sceptre and ondisc R-packages are installed. The script
## was tested with Java 23.0.2 installed.

##########################
# REQUIRED INPUT ARGUMENTS
##########################
data_directory="/g/steinmetz/project/otar_2063/promoter_screen/processing/combine_data/results/combined_sceptre_filtered"
# sceptre object
sceptre_object_fp=$data_directory"/sceptre_object.rds"
# response ODM
response_odm_fp=$data_directory"/gene.odm"
# grna ODM
grna_odm_fp=$data_directory"/grna.odm"
# object containing model formula
formula_object="./full_covars_formula.rds"

###################
# OUTPUT DIRECTORY:
##################
output_directory="./sceptre_outputs"

#################
# Invoke pipeline
#################
nextflow run timothy-barry/sceptre-pipeline -r main \
  --sceptre_object_fp $sceptre_object_fp \
  --response_odm_fp $response_odm_fp \
  --grna_odm_fp $grna_odm_fp \
  --output_directory $output_directory \
  --formula_object $formula_object \
  --grna_assignment_formula  $formula_object \
  --grna_pod_size 250 \
  --assign_grnas_time_per_grna 180s \
  --assign_grnas_memory 24GB \
  --pipeline_stop assign_grnas \
  -resume
