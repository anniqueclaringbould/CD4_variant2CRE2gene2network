
# load filtering parameters
configfile: "config/alternative_filters.yml"

# alternative filtering
rule alternative_filter:
  input:
    etps = "results/results_df_with_promoterC_annotated_allFeatures.csv",
    gene_annot = "resources/genes.gtf"
  output:
    pairs = "results/alternative_filters/{filter_scenario}/all_pairs.tsv.gz",
    positives = "results/alternative_filters/{filter_scenario}/positives.tsv.gz",
    negatives = "results/alternative_filters/{filter_scenario}/negatives.tsv.gz"
  params:
    filters = lambda wildcards: config["alternative_filters"][wildcards.filter_scenario]
  conda: "../envs/otar_cre_figures.yml"
  resources:
    mem_mb = 32000
  script:
    "../scripts/alternative_filter.R"

# compare alternative filters
rule plot_filter_scenarios:
  input:
    positives = expand("results/alternative_filters/{filter_scenario}/positives.tsv.gz",
      filter_scenario = config["alternative_filters"]),
    negatives = expand("results/alternative_filters/{filter_scenario}/negatives.tsv.gz",
      filter_scenario = config["alternative_filters"]),
    scenarios = "config/filter_scenarios.tsv"
  output: "results/plots/alternative_cre_screen_filters.pdf"
  conda: "../envs/plot_filter_scenarios.yml"
  script:
    "../scripts/plot_filter_scenarios.R"
