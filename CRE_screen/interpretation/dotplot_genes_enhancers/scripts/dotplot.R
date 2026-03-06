####################################################################
## Script to make dot plot of enhancers and target genes by panel ##
####################################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(openxlsx)
library(pals)

#directories
root <- "/g/steinmetz/project/otar/"
res_dir <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM/results/"
outdir <- paste0(root, "interpretation/dotplot_genes_enhancers/output/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "select_genes_full_screen/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read results from all 16 panels
combined_results <- list()

for (i in 1:16) {
  file_path <- paste0(res_dir, "Panel", i, "/discovery_results.tsv")
  
  if (file.exists(file_path)) {
    res <- fread(file_path)
    res <- res %>% mutate(panel = i)  # Add a column with the panel number
    combined_results[[i]] <- res
  } else {
    message(paste("File not found for panel", i, ":", file_path))
  }
}

#combine all results into one data frame
res <- do.call(rbind, combined_results)

#read gene definition file
gtf <- readGFF("/g/steinmetz/brausche/genomes/alias/hg38/gencode_gtf/default/hg38.gtf")

#read chromosome size file
chromosome <- fread(paste0(root, "interpretation/dotplot_genes_enhancers/data/chr_sizes.txt"))

#### ADJUST ####

#keep gene name to ENSG conversion and gene location
gtf <- gtf %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name, seqid, start, end) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  distinct() 

#add ENSG IDs and gene locations to res, and change control enhancer/promoter names to their locations
res <- res %>%
  left_join(gtf, by = c("response_id" = "gene_name")) %>% #add ENSG IDs and positions
  mutate(grna_target = case_when(grna_target == "CD28" ~ "chr2:203706639-203706639", #add promoter controls as gene names (as in results)
                                 grna_target == "CTLA4" ~ "chr2:203867771-203867771",
                                 grna_target == "IL2RA" ~ "chr10:6062367-6062367",
                                 grna_target == "CD81" ~ "chr11:2377310-2377310",
                                 grna_target == "CD4" ~ "chr12:6789528-6789528",
                                 grna_target == "SP140" ~ "chr2:230225130-230226457", #add pilot controls as gene names (as in results)
                                 grna_target == "CAB39" ~ "chr2:230658624-230659752",
                                 grna_target == "PIM1" ~ "chr6:37050015-37051033",
                                 grna_target == "CCND2" ~ "chr12:4153422-4154015",
                                 grna_target == "JUND" ~ "chr19:18291227-18293782",
                                 .default = grna_target)) %>%
  separate(grna_target, remove = F, into = c("enhancer_chr", "enhancer_start", "enhancer_end")) %>% #save enhancer position
  mutate(link_id = paste0(grna_target, "_", gene_id)) #save link IDs

#filter results to keep only links that pass QC
res2 <- res %>%
  filter(pass_qc == T) %>% #keep only results that pass QC
  select(link_id, grna_target, enhancer_chr, enhancer_start,
         gene_name = response_id, gene_id, gene_chr = seqid, gene_start = start,
         p_value, log_2_fold_change, significant, panel) %>%
  group_by(link_id) %>%
  slice(which.min(p_value)) %>% #for each link, keep only one row with results based on minimal p-value 
  ungroup()

#Adjust links data
links2 <- links %>%
  mutate(SNP_trait = case_when(SNP_trait == "" ~ source, #add source to SNP trait for non GWAS links
                               .default = SNP_trait),
         peak_type = case_when(peak_type == "" ~ "control", #add peak type to peaks that were selected as controls
                               .default = peak_type)) %>%
  filter(SNP_trait != "nonGWAS_control_promoter") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_control_pilot") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_control_gsk") %>% #these are in there double
  filter(source != "aim1") %>% #remove the genes that were added for aim1, since they have no enhancer connected
  mutate(link_id = paste0(peak, "_", gene_id))

#Keep the link IDs of prioritised links
prioritised_links <- links2 %>% 
  separate(peak, remove = F, into = c("enhancer_chr", "enhancer_start", "enhancer_end")) %>% #save enhancer position
  left_join(gtf, by = c("gene_id", "gene_name")) %>% #add chromosome positions of the genes
  select(link_id, grna_target = peak, enhancer_chr, enhancer_start,
         gene_name, gene_id, gene_chr = seqid, gene_start = start, gene_end = end) %>% 
  distinct()

#Keep link IDs of links that failed QC
failed_QC_links <- res %>% 
  filter(pass_qc == F) %>% #keep only results that fail QC
  select(grna_target, gene_id, link_id) %>% 
  distinct()

#Add absolute coordinates to chromosome sizes
chromosome <- chromosome %>% 
  mutate(abs_coord = cumsum(as.numeric(Total_length_bp)), #make cumulative sum minus the first one since it should be the start point of each chr
         tick_mark = abs_coord - (Total_length_bp/2)) %>%
  filter(Chromosome != "Y")

max_chr <- max(chromosome$abs_coord)

res3 <- res2 %>%
  full_join(prioritised_links, by = c("link_id", "grna_target", "enhancer_chr", "enhancer_start", "gene_name", "gene_id", "gene_chr", "gene_start")) %>% #add links that may not have passed QC or that were not tested but they were prioritised
  mutate(link_type = case_when(link_id %in% prioritised_links$link_id & significant == T ~ "prioritised, significant",
                               link_id %in% prioritised_links$link_id & significant == F ~ "prioritised, not significant",
                               !link_id %in% prioritised_links$link_id & significant == T ~ "not prioritised, significant",
                               !link_id %in% prioritised_links$link_id & significant == F ~ "not prioritised, not significant",
                               link_id %in% prioritised_links$link_id & link_id %in% failed_QC_links$link_id ~ "prioritised, failed QC",
                               link_id %in% prioritised_links$link_id & !link_id %in% res$link_id ~ "prioritised, not tested",
                               .default = "missing")) %>%
  mutate(enhancer_chr = str_replace(enhancer_chr, "chr", ""),
         gene_chr = str_replace(gene_chr, "chr", ""),
         enhancer_start = as.integer(enhancer_start),
         panel = as.factor(panel)) %>%
  select(link_id, enhancer_chr, enhancer_start, gene_chr, gene_start, significant, link_type, p_value, panel)

table(res3$link_type, useNA = 'always')

# Calculate absolute chromosome coordinates for enhancers
res3$enhancer_abspos <- 0
res3[res3$enhancer_chr == 1, ]$enhancer_abspos <- res3[res3$enhancer_chr == 1, ]$enhancer_start
for (i in 2:22){
  res3[res3$enhancer_chr == i, ]$enhancer_abspos <- res3[res3$enhancer_chr == i, ]$enhancer_start + chromosome[chromosome$Chromosome == i - 1, ]$abs_coord
}
res3[res3$enhancer_chr == "X", ]$enhancer_abspos <- res3[res3$enhancer_chr == "X", ]$enhancer_start + chromosome[chromosome$Chromosome == 22, ]$abs_coord

# Calculate absolute chromosome coordinates for genes
res3$gene_abspos <- 0
res3[res3$gene_chr == 1, ]$gene_abspos <- res3[res3$gene_chr == 1, ]$gene_start
res3$gene_abspos <- as.numeric(res3$gene_abspos)
for (i in c(2:22)){
  res3[res3$gene_chr == i, ]$gene_abspos <- res3[res3$gene_chr == i, ]$gene_start + chromosome[chromosome$Chromosome == i - 1, ]$abs_coord
}
res3[res3$gene_chr == "X", ]$gene_abspos <- res3[res3$gene_chr == "X", ]$gene_start + chromosome[chromosome$Chromosome == 22, ]$abs_coord

#Add variables for colouring and alpha values
res4 <- res3 %>%
  mutate(colour_var = factor(case_when(link_type == "prioritised, significant" ~ paste0("Panel ", panel, " prioritised hit"),
                                       link_type == "not prioritised, significant" ~ "unexpected hit",
                                       .default = "not significant"),
                             levels = c(paste0("Panel ", 1:16, " prioritised hit"), "unexpected hit", "not significant")),
         alpha_var = case_when(colour_var == "not significant" ~ "not significant",
                               colour_var == "unexpected hit" ~ "unexpected hit",
                               .default = "significant"),
         p_value = case_when(is.na(p_value) ~ 1, #add p-value of 1 for links that were not in the results but were prioritised
                             .default = p_value)) 

##### PLOT #####

#frame / background raster
p <- ggplot(chromosome, aes(x = abs_coord, y = abs_coord)) + 
  theme_bw() + 
  geom_vline(xintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + 
  geom_hline(yintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + 
  scale_x_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$Chromosome, sep = ''), limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$Chromosome, sep = ''), limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  xlab('\nEnhancer position (hg38)') + 
  ylab('Gene start (hg38)\n')

#use different alpha values
my_alphas <- c("not significant" = 0.1,
               "unexpected hit" = 0.3,
               "significant" = 0.5)

p + geom_point(data = as.data.frame(res4[res4$colour_var == "not significant",]), #add not significant results as grey dots in the background
               aes(x = enhancer_abspos,
                   y = gene_abspos,
                   size = -log10(p_value),
                   colour = colour_var,
                   alpha = alpha_var)) +
  geom_point(data = as.data.frame(res4[res4$colour_var != "not significant",]), #overlay significant results
             aes(x = enhancer_abspos,
                 y = gene_abspos,
                 size = -log10(p_value),
                 colour = colour_var,
                 alpha = alpha_var)) +
  scale_colour_manual(values=c("grey90", unname(stepped())[1:17])) +
  scale_alpha_manual(values = my_alphas) +
  scale_size_continuous(breaks = c(1, 10, 100, 300), range = c(0.75 * 1, 0.75 * 5),
                        guide = guide_legend(title = expression(paste(-log[10]("P-value"))),
                                             override.aes = list(alpha = 1))) +
  guides(alpha = "none",
         colour = guide_legend(title = ""))

ggsave(paste0(outdir, 'dot_plot_enhancers_genes_hits_highlighted_per_panel', Sys.Date(), '.pdf'), height = 9, width = 9 * 1.25, dpi = 400)


