
# download chromatin assay files
rule download_chromatin_assay:
  output: "resources/chromatin_assays/{file}.bigWig"
  conda: "../envs/otar_cre_figures.yml"
  shell:
    "wget -O {output} https://www.encodeproject.org/files/{wildcards.file}/@@download/{wildcards.file}.bigWig"

## Chromatin activity assays -----------------------------------------------------------------------
    
# extract enhancer elements in .bed format
rule extract_elements:
  input: "results/results_df_with_promoterC_annotated.csv"
  output: temp("results/enhancers.bed")
  params:
    elements = "enhancers"  
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/extract_elements.R"
    
# extract promoter regions in .bed format
rule extract_promoters:
  input: "results/results_df_with_promoterC_annotated.csv"
  output: temp("results/promoters.bed")
  params:
    elements = "promoters",
    extend_tss = 250
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/extract_elements.R"

# compute chromatin assay signal at elements 
rule count_chromatin_assay_bigwig:
  input:
    assay = "resources/chromatin_assays/{file}.bigWig",
    elements = "results/{elements}.bed",
  output: "results/chromatin_assays/{elements}.{file}.tsv"
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/count_chromatin_assay_bigwig.R"

# count and combine signal for all assays and metrics, and add to E-G results table
rule combine_chromatin_assays:
  input:
    enh_assays = expand("results/chromatin_assays/enhancers.{file}.tsv",
      file = config["chromatin_assays"].values()),
    prom_assays  = expand("results/chromatin_assays/promoters.{file}.tsv",
      file = ["ENCFF546IYU", "ENCFF600IHC"]),
    elements = "results/results_df_with_promoterC_annotated.csv"
  output: temp("results/results_df_with_promoterC_annotated_chromFeatures.csv")
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/combine_chromatin_assays.R"
  
# add additional functional annotations, like expression level, eRNA counts or disease information
rule add_additional_functional_data:
  input:
    elements = "results/results_df_with_promoterC_annotated_chromFeatures.csv",
    eRNA_counts = config["eRNA_counts"],
    gene_tpm = "resources/gene_tpm.csv.gz",
    gene_constraints = "resources/gene_constraints.csv.gz",
    ubiq_expression = "resources/Ubiquitous_expression.txt.gz",
    gene_disease_info = config["gene_disease_info"]
  output: temp("results/results_df_with_promoterC_annotated_chromFeaturesPlusOthers.csv")
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/add_additional_functional_data.R"
    
# calculate Hi-C interaction frequencies between E-G pairs on the same chromsomes
rule add_hic_interaction_frequency:
  input:
    elements = "results/results_df_with_promoterC_annotated_chromFeaturesPlusOthers.csv",
    hic = config["hic"]["file"]
  output: "results/results_df_with_promoterC_annotated_allFeatures.csv"
  params:
    hic_res = config["hic"]["resolution"]
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/compute_hic_contacts.py"
