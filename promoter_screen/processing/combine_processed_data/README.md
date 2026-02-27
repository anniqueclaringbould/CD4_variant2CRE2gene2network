# Combine processed data
This snakemake workflow combines processed sequencing data from different samples into one Seurat
object using memory-efficient on disk storage implemented in the BPCells R packages. It also creates
the SCEPTRE and AnnData objects required for performing guide assignments and differential
expression tests.
