# Locus network plots

This code was used to create the networks for the two example loci shown in Figure 5. To run these,
outputs from both the `CRE_screen/interpretation/chromatin_analyses` and 
`promoter_screen/interpretation/DEG_testing` workflows are needed.

The `1_get_promoter_screen_hits.R` is used to extract significant hits from the promoter screen
differential expression testing results. `2_make_locus_networks_for_manuscript.R` then creates locus
network plots using functions defined in `locus_network_functions.R`.
