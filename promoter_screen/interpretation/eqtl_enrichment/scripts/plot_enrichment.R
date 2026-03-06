####################################
## Script to plot eQTL enrichment ##
## cis-trans links from eQTLGen p2 #
####################################

#### SETUP ####
#libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rtracklayer)

#directories
root    <- "/g/steinmetz/project/otar_2063/"
eqtldir <- paste0(root, "promoter_screen/interpretation/eqtl_enrichment/")
outdir  <- paste0(eqtldir, "output/")

#### READ DATA ####
pairs <- fread(paste0(outdir, "promoter_screen_with_overlapping_cis-trans-eQTLs_pseudobulk_deseq.tsv"))
tab <- fread(paste0(outdir, "contingency_table_gene-gene_pairs_pseudobulk_deseq.tsv"))

#### ADJUST ####
tab_pct <- tab %>%
  pivot_wider(names_from = eqtl_annotated, names_prefix = "eqtl_annotated_", values_from = Freq) %>%
  mutate(pct_overlap = eqtl_annotated_TRUE/eqtl_annotated_FALSE,
         sig_label = case_when(significant == TRUE ~ "significant\n gene pairs",
                               significant == FALSE ~ "not significant\n gene pairs"),
         sig_label = factor(sig_label, levels = c("not significant\n gene pairs", "significant\n gene pairs")),
         label = paste0(eqtl_annotated_TRUE, " / ", eqtl_annotated_FALSE))

#### STATS ####
fisher <- tab %>%
  pivot_wider(names_from = eqtl_annotated, names_prefix = "eqtl_annotated_", values_from = Freq) %>%
  column_to_rownames(var = "significant") %>%
  fisher.test(.)

label_text <- paste0("OR = ", round(fisher$estimate, 2), ", p = ", signif(fisher$p.value, 3))

#### PLOT ####
ggplot(tab_pct, aes(x = sig_label, y = pct_overlap)) +
  geom_bar(aes(fill = significant), col = "black",  stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.001)) +
  scale_fill_manual(values = c("TRUE" = "#0263a3", "FALSE" = "grey60")) +
  labs(x = NULL,
    y = "\n% overlap with eQTL-derived gene pairs") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_bracket(xmin = 1, xmax = 2,
               size = 0.8,
               label = label_text,
               y.position = 1.05 * max(tab_pct$pct_overlap),
               vjust = 4.8,
               hjust = 1)

ggsave(paste0(outdir, "eQTL_enrichment_cis-trans_gene_pairs.pdf"), height = 2, width = 5)
 