## Other enhancer screen related analyses

# add enhancer screen guide assignments and sceptre objects to config
configfile: "config/enhancer_screen_sceptre_files.yml"

# make enhancer screen volcano plot and color pairs by prioritization strategy
rule prioritization_volcano_plot:
  input:
    enh_results = "results/results_df_with_promoterC_annotated_allFeatures.csv",
    pr_classes = config["eg_prioritization"]
  output: "results/plots/prioritization_volcano_plot.pdf"
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/prioritization_volcano_plot.R"

# make relative expression levels per TAP-seq target gene panel
rule tapseq_panel_summary_plots:
  input:
    grna_targets = "results/guide_assignments/grna_targets.tsv.gz",
    guide_assignments = "results/guide_assignments/guide_assignments.tsv.gz",
    gene_panels = config["tapseq_gene_panels"],
    gene_tpm = "resources/gene_tpm.csv.gz",
    enh_results = "results/results_df_with_promoterC_annotated_allFeatures.csv"
  output: "results/plots/tapseq_panel_summary_plots.pdf"
  params:
    filt_assigned_guides = False,
    filt_detected_genes = False
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/tapseq_panel_summary_plots.R"
    
# make sceptre QQ plots for all CRE screen panels
rule make_cre_screen_qq_plots:
  input: config["finished_sceptre_objects"].values()
  output:
    qq_plots_tail_2by8 = "results/plots/cre_screen_sceptre_qq_plots_tail_2by8.pdf",
    qq_plots_tail_4by4 = "results/plots/cre_screen_sceptre_qq_plots_tail_4by4.pdf"
  conda: "../envs/otar_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/make_cre_screen_qq_plots.R"    

# process enhancer screen guide assignments for one panel
rule process_guide_assignments:
  input:
    guide_assignments = lambda wildcards: config["guide_assignments"][wildcards.panel],
    sceptre_object = lambda wildcards: config["sceptre_objects"][wildcards.panel]
  output:
    guide_assignments = temp("results/guide_assignments/guide_assignments_{panel}.tsv.gz"),
    summary_stats = temp("results/guide_assignments/summary_stats_{panel}.tsv.gz"),
    grna_targets = temp("results/guide_assignments/grna_targets_{panel}.tsv.gz")
  conda: "sceptre"  # environment needs to be built separately
  resources:
    mem = "16G"
  script:
    "../scripts/miscellaneous_analyses/process_guide_assignments.R"
    
# combine enhancer screen guide assignments into one file
rule combine_guide_assignments:
  input:
    guide_assignments = expand("results/guide_assignments/guide_assignments_{panel}.tsv.gz",
      panel = config["guide_assignments"]),
    summary_stats = expand("results/guide_assignments/summary_stats_{panel}.tsv.gz",
      panel = config["guide_assignments"]),
    grna_targets = expand("results/guide_assignments/grna_targets_{panel}.tsv.gz",
      panel = config["guide_assignments"])
  output: 
    guide_assignments = "results/guide_assignments/guide_assignments.tsv.gz",
    summary_stats = "results/guide_assignments/summary_stats.tsv.gz",
    grna_targets = "results/guide_assignments/grna_targets.tsv.gz"
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/combine_guide_assignments.R"

# calculate miscellaneous enhancer screen numbers
rule calculate_misc_numbers:
  input:
    guide_assignments = "results/guide_assignments/guide_assignments.tsv.gz",
    summary_stats = "results/guide_assignments/summary_stats.tsv.gz"
  output: "results/misc_numbers.txt"
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/calculate_misc_numbers.R"
    
# make results table to share with OTAR
rule make_sharable_results_table:
  input: "results/results_df_with_promoterC_annotated_allFeatures.csv"
  output: "results/share/cre_perturbation_screen_results.csv"
  conda: "../envs/otar_cre_figures.yml"
  script:
    "../scripts/miscellaneous_analyses/make_sharable_results_table.R"
