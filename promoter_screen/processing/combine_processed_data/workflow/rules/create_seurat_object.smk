## Create combined Seurat object using BPCells' on-disk functionality 

# combine metric summary files from all samples
rule combine_metric_summary_files:
  input: config["metric_files"]
  output: "results/combined_metric_summaries.rds"
  params:
    metrics = config["qc_metrics"]
  conda: "../envs/combine_data.yml"
  script:
    "../scripts/combine_metric_summary_files.R"

# combine seurat files for one donor to help mitigate memory usage
rule combine_seurat_donor:
  input: lambda wildcards: config["seurat_files"][wildcards.donor]
  output: 
    seu = "results/per_donor_and_chip/{donor}_combined_seurat.rds",
    counts_dir = directory("results/per_donor_and_chip/{donor}_combined_counts")
  log: "results/logs/combine_seurat_{donor}.log"
  conda: "../envs/combine_data.yml"
  resources:
    mem = "256G",
    runtime = "3h",
    partition = "bigmem"
  script:
    "../scripts/combine_seurat_donor.R"

# combine seurat objects from all donors and convert to on-disk object
rule combine_seurat_all_donors:
  input:
    expand("results/per_donor_and_chip/{donor}_combined_seurat.rds", donor = config["seurat_files"])
  output: 
    seu = "results/combined_seurat.rds",
    counts_dir = directory("results/combined_counts")
  log: "results/logs/combine_seurat_all_donors.log"
  conda: "../envs/combine_data.yml"  
  resources:
    mem = "24G",
    runtime = "8h"
  script:
    "../scripts/combine_seurat_all_donors.R"
