
# plot chromatin assays per broad ETP class
rule chromatin_features_etp_classes:
  input: "results/results_df_with_promoterC_annotated_allFeatures.csv"
  output: "results/chromatin_features_etp_classes.html"
  params:
    hic_res = config["hic"]["resolution"]
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/chromatin_features_etp_classes.Rmd"
    
# analyze chromatin landscape of enhancer- vs. silencer-like regualtory elements and interactions
rule enhancer_vs_silencer_analysis:
  input: 
    annot_pairs = "results/results_df_with_promoterC_annotated_allFeatures.csv",
    all_pairs = "results/results_df_with_promoterC_annotated.csv",
    gene_tpm = "resources/gene_tpm.csv.gz",
    annot = "resources/hg38.gtf",
    ray2025_results = "resources/other_datasets/ray_et_al_2025_table_S3.tsv.gz",
    gasperini2019_results = "resources/other_datasets/ENCODE_Gasperini2019_0.13gStd_MAST_perCRE_GRCh38.tsv.gz"
  output: "results/enhancer_vs_silencer_analysis.html"
  conda: "../envs/otar_cre_figures.yml"
  resources:
    mem_mb = 32000
  script:
    "../scripts/enhancer_vs_silencer_analysis.Rmd"
