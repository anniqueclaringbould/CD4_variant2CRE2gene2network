
# plot chromatin assays per broad ETP class
rule chromatin_features_etp_classes:
  input: "results/results_df_with_promoterC_annotated_allFeatures.csv"
  output: "results/chromatin_features_etp_classes.html"
  params:
    hic_res = config["hic"]["resolution"]
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/chromatin_features_etp_classes.Rmd"
