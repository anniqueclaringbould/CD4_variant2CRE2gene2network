## create genes reference for downstream analyses

# save.image("genes.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(rtracklayer)
})

# load genome annotations and TAP-seq target genes
annot <- import(snakemake@input$genome)
target_genes <- import(snakemake@input$target_genes)

# only retain genes on autosomes and chrX
chrs <- paste0("chr", c(1:22, "X"))
annot <- annot[seqnames(annot) %in% chrs]

# extract exons of protein-coding and lncRNA genes
gene_types <- c("protein_coding", "lincRNA")
exons <- annot[annot$type == "exon" &
                annot$gene_type %in% gene_types &
                annot$transcript_type %in% gene_types]

# strip version from gene_ids
exons$gene_id <- sub("\\..*", "", exons$gene_id)

# get annotations on genes not part of the target genes
other_genes <- setdiff(unique(exons$gene_id), unique(target_genes$gene_id))
nontarget_genes <- exons[exons$gene_id %in% other_genes]

# combine with target genes annotations to get genome-wide annotations of gene exons
genes <- sort(c(target_genes, nontarget_genes))

# save to output file
export(genes, con = snakemake@output[[1]])
