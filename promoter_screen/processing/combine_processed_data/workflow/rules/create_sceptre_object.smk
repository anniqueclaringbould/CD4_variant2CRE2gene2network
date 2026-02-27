## Create sceptre object from QC filtered seurat object

# Create all required input files ------------------------------------------------------------------

# create QC filtered counts matrices
rule write_filt_mtx_donor:
  input:
    seu = config["filtered_seurat"],
    annot = config["annot"],
    rm_feat = "resources/guides_to_remove.txt"
  output:
    mtx = "results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/matrix.mtx.gz",
    barcodes = "results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/barcodes.tsv.gz",
    features = "results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/features.tsv.gz"
  log: "results/logs/write_filt_mtx_{donor}.log"
  params:
    grna_pattern = config["grna_pattern"]
  conda: "../envs/combine_data.yml"
  resources:
    mem = "64G",
    runtime = "4h",
    partition = "bigmem"
  script:
    "../scripts/write_mtx_donor.R"
    
# create additional covariate matrix for filtered sceptre object
rule create_sceptre_covariates:
  input: 
    seu = config["filtered_seurat"],
    add_covars = "resources/cell_state_covariates.rds"
  output: "results/combined_sceptre_filtered/cell_covariates.rds"
  params:
    seurat_covariates = {"percent_rb": "percent.rb", "percent_mt": "percent.mt"},
    seurat_pcs = [1, 2, 3, 4, 5, 6]
  conda: "../envs/combine_data.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/create_sceptre_covariates.R"
    
# create gRNA targets table based on gRNA ids in UMI count features
rule create_sceptre_grna_targets:
  input:
    expand("results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/features.tsv.gz",
           donor = config["seurat_files"])
  output: "results/combined_sceptre_filtered/grna_targets.tsv.gz"
  conda: "../envs/combine_data.yml"
  script:
    "../scripts/create_sceptre_grna_targets.R"
    
  # create a list of self promoter targeting perturbations to use as positive control pairs    
rule create_sceptre_positive_control_pairs:
  input:
    grna_targets = "results/combined_sceptre_filtered/grna_targets.tsv.gz",
    features = expand("results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/features.tsv.gz",
                      donor = config["seurat_files"])
  output: "results/combined_sceptre_filtered/positive_control_pairs.tsv.gz"
  conda: "../envs/combine_data.yml"
  script:
    "../scripts/create_sceptre_positive_control_pairs.R"  

# Create sceptre objects ---------------------------------------------------------------------------

# create filtered on-disc sceptre object for data from all donors
rule create_sceptre_object:
  input:
    mtx = expand("results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/matrix.mtx.gz",
                  donor = config["seurat_files"]),
    barcodes = expand("results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/barcodes.tsv.gz",
                      donor = config["seurat_files"]),
    features = expand("results/per_donor_and_chip/{donor}_combined_counts_filt_mtx/features.tsv.gz",
                      donor = config["seurat_files"]),
    grna_targets = "results/combined_sceptre_filtered/grna_targets.tsv.gz",
    covariates = "results/combined_sceptre_filtered/cell_covariates.rds",
    pos_ctrl_pairs = "results/combined_sceptre_filtered/positive_control_pairs.tsv.gz"
  output:
    "results/combined_sceptre_filtered/sceptre_object.rds",
    "results/combined_sceptre_filtered/gene.odm",
    "results/combined_sceptre_filtered/grna.odm"
  log: "results/logs/create_filt_sceptre_object.log"
  params:
    moi = config["sceptre"]["moi"],
    side = config["sceptre"]["side"],
    resampling_mechanism = config["sceptre"]["resampling_mechanism"],
    formula = config["sceptre"]["formula"],
    on_disc = True
  conda: "sceptre"
  resources:
    mem = "48G",
    runtime = "6h",
  script:
    "../scripts/create_sceptre_object.R"
    
# create filtered on-disc sceptre object for data from one donor only
rule create_sceptre_object_donor:
  input:
    mtx = lambda wildcards: expand("results/per_donor_and_chip/{sample}_combined_counts_filt_mtx/matrix.mtx.gz",
                                    sample = config["samples_per_donor"][wildcards.donor]),
    barcodes = lambda wildcards: expand("results/per_donor_and_chip/{sample}_combined_counts_filt_mtx/barcodes.tsv.gz",
                                        sample = config["samples_per_donor"][wildcards.donor]),
    features = lambda wildcards: expand("results/per_donor_and_chip/{sample}_combined_counts_filt_mtx/features.tsv.gz",
                                        sample = config["samples_per_donor"][wildcards.donor]),
    grna_targets = "results/combined_sceptre_filtered/grna_targets.tsv.gz",
    covariates = "results/combined_sceptre_filtered/cell_covariates.rds",
    pos_ctrl_pairs = "results/combined_sceptre_filtered/positive_control_pairs.tsv.gz"
  output:
    "results/combined_sceptre_filtered_per_donor/{donor}/sceptre_object.rds",
    "results/combined_sceptre_filtered_per_donor/{donor}/gene.odm",
    "results/combined_sceptre_filtered_per_donor/{donor}/grna.odm"
  log: "results/logs/create_filt_sceptre_object_donor_{donor}.log"
  params:
    moi = config["sceptre"]["moi"],
    side = config["sceptre"]["side"],
    resampling_mechanism = config["sceptre"]["resampling_mechanism"],
    formula = config["sceptre"]["formula"],
    on_disc = True
  conda: "sceptre"
  resources:
    mem = "48G",
    runtime = "6h",
  script:
    "../scripts/create_sceptre_object.R"
