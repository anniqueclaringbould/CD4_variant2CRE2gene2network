# Promoter screen sequencing data processing

This snakemake workflow aligns promoter screen sequencing data to the a customized genome containing
the human genome and used gRNA sequences, and extracts gene and gRNA UMI counts for each cell using
the BD Rhapsody Sequence Analysis Pipeline.

Running this workflow requires the
[BD Rhapsody Sequence Analysis Pipeline](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline) is
installed and available via the systems path. The workflow was used with version 2.2 of the BD
Rhapsody pipeline.

A generalizable template of this workflow for new datasets is available here: https://github.com/argschwind/ReptilianRhapsody
