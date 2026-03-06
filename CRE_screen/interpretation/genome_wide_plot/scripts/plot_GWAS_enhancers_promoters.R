###############################################################################
## Script to make Manhattan plot with GWAS signals, enhancers, and promoters ##
###############################################################################

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
library(pheatmap)
library(ieugwasr)
library(ggrastr)
library(patchwork)
library(RColorBrewer)

#directories
root <- "/g/steinmetz/project/otar/"
outdir <- paste0(root, "interpretation/genome_wide_plot/output/")
enhdir <- paste0(root, "select_enhancer_peaks/")
genedir <- paste0(root, "select_genes_full_screen/")
promdir <- paste0(root, "promoter_screen/output/")
  
#### READ DATA ####

#GWAS SNPs we included in the enhancer screen
snps <- fread(paste0(enhdir, "Final_SNPs_large_screen_2024-03-01.bed"))

#GWA studies
gwas <- fread(paste0(root, "data/gwas_selection/filtered_studies_OpenGWAS_2023-02-10.tsv"))

#Enhancers we included in the enhancer screen
enhs <- fread(paste0(enhdir, "Final_peaks_large_screen_2024-03-01.tsv"))

#Enhancer activities
enh_act <- fread(paste0(enhdir, "normalised_peaks_across_conditions/ATAC_counts_normToMax_quantileNorm_euclideanNorm_extended_1000_filtered.txt"))

#Genes we included in the enhancer screen
genes <- fread(paste0(root, "select_genes_full_screen/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_unique_genes_2024-03-01.tsv"))

#Genes we included in the promoter screen
proms <- fread(paste0(promdir, "gene_list_filtered_protein_coding_TPM10_TFs_enhancer_targets.txt"))

#All genes in T cells
genes_exp <- fread(paste0(promdir, "gene_list_unfiltered.txt"))

#TFs we included in the promoter screen
tf <- fread("/g/zaugg/zaugg_shared/annotations/TFBS/hg38/PWMScan_HOCOMOCOv12/HOCOMOCOv12INVIVO_SCENICplus_mapping.tsv")

# #all gwas snps that we initially selected
# gwas_snps_files <- list.files(path = paste0(root, "data/gwas_selection/"), pattern = "snplist.tsv", recursive = T, full.names = T)
# fread_snps <- function(file){
#   disease <- basename(dirname(file))  # Extract folder name (disease)
#   dm <- fread(file, header = F) %>%
#     rename(rsid = V1) %>%
#     mutate(disease = disease)
# }
# gwas_snps <- do.call(rbind, lapply(gwas_snps_files, fread_snps))

#functions
select = dplyr::select
mutate = dplyr::mutate

#### ADJUST ####

#change disease names and keep only conversion info
gwas_conv <- gwas %>%
  mutate(trait = str_replace_all(trait, " ", "_")) %>%
  mutate(trait = str_replace_all(trait, "'", "")) %>%
  select(id, trait)

#add trait groups
disease_groups <- c("Allergic_disease" = "Skin and allergic",
                    "Asthma" = "Skin and allergic",
                    "Eczema" = "Skin and allergic",
                    "Psoriasis" = "Skin and allergic",
                    "Celiac_disease" = "Gastrointestinal",
                    "Crohns_disease" = "Gastrointestinal",
                    "Inflammatory_bowel_disease" = "Gastrointestinal",
                    "Ulcerative_colitis" = "Gastrointestinal",
                    "Ankylosing_spondylitis" = "Rheumatic",
                    "Gout" = "Rheumatic",
                    "Rheumatoid_arthritis" = "Rheumatic",
                    "Systemic_lupus_erythematosus" = "Rheumatic",
                    "Multiple_sclerosis" = "Rheumatic",
                    "Primary_biliary_cholangitis" = "Hepatic",
                    "Primary_biliary_cirrhosis" = "Hepatic",
                    "Primary_sclerosing_cholangitis" = "Hepatic",
                    "Type_1_diabetes" = "Type_1_diabetes")

#rename SNPs we ended up including and add study ID
names(snps) <- c("chr", "start", "end", "rsid", "trait")
snps <- snps %>%
  left_join(gwas_conv, by = "trait") %>%
  mutate(id = case_when(is.na(id) ~ "nonGWAS",
                        .default = id),
         disease_group = disease_groups[trait])

#make a list of studies 
studies <- unique(snps[snps$id != "nonGWAS",]$id)

#save token to talk to OpenGWAS database
#opengwas_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJhbm5pcXVlY2xhcmluZ2JvdWxkQGdtYWlsLmNvbSIsImlhdCI6MTcyNTk3MjU2MCwiZXhwIjoxNzI3MTgyMTYwfQ.qk7a87O3K4N5NnhJVqyk8qcxRXJHIe8kk8zVO6Ik2q2C2DGCw-N-8MGVNZb5PKSUaqlRFVz6UkhUn0qwGLoAyxP5jNWdl3Urhh2CDEqS0iJnKmOouxo1gsq_8DrGwhfw4O0H9USV1fJGTXGDpD01RWL1yQh1hevI9mgDbiyZzPkV16hoXEboGpl1Jr2tD_xNrUZt9aFYt2H56uxEz-_ZlnqgjzdVdfKd_WAT9VWBQCe_fgy2_ZWZf6hO6akHgjjcy92mFXlZ7KOA4pPAR4tn7Tui4s_YPIH6HuE6PUhF6EY35psk_c3lgI0eQrYHC2yetUKd4KGDFifcNMmYl3ZrhQ"
opengwas_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJhbm5pcXVlY2xhcmluZ2JvdWxkQGdtYWlsLmNvbSIsImlhdCI6MTcyNzQzOTc3NiwiZXhwIjoxNzI4NjQ5Mzc2fQ.WB7FR3EVh5IhLoQpaXx1q_3_e9Zo_DoRkAlNKmGCDOZA_dxv9LHb18_95F-RBSpm5D3FfyHJhoLvoMh6-23AE18ZwiWWV4OKGogjmQL_zkdkw9NmUBjG87w0FRoBXjWlX1ChzJ9Npe6YFDVfqJAm7oSNYZTlqeXaXy_S45fgtzz-rptX7aH4Y79pXdM_gA_Q5h58BVBEpTmYV_ZPXN92IfvEgY3djj-dQtmIcG2UA3cw_Ot2v8lVrH21GYFeMB31C7wvasJHQSBprWmOb425QfEIGFv8B6Mc8mU-Vy5d3XJZ0b0s44c_Hl8WQO8ozYWx9OZs8vAUSzIg_VtnLnGkWQ"

#Get all GWAS results (with p-values) for the diseases we ended up including
#Specifically query the SNPs we ended up including
get_hits <- function(study){
  cat("Now submitting", study, "\n")
  
  #find almost all variants per GWAS for the Manhattan plot (remove p < 0.05 to keep the size manageable)
  all <- tophits(id = study,
               opengwas_jwt = opengwas_jwt,
               clump = 0, #no clumping
               pval = 0.05)
  
  #specifically look for the SNPs we included through CHEERS
  snps_disease <- snps %>%
    filter(id == !!study) %>%
    pull(rsid)
    
  selected <- associations(variants = snps_disease, 
                    id = study,
                    opengwas_jwt = opengwas_jwt)
  
  #keep the same names as the 'all' data frame
  selected <- selected %>%
    select(se, p, chr, position, beta, n, id, rsid, ea, nea, eaf, trait) 
  
  cat(length(unique(selected$rsid)), "of", length(unique(snps_disease)), "query SNPs found\n")
  
  d <- rbind(all, selected) %>%
    distinct()
  
  cat(length(unique(d$rsid)), "total SNPs found\n")
  
  d <- d %>%
    mutate(trait = trait %>%
             str_replace_all(" \\(.*", "") %>%
             str_replace_all(" ", "_") %>%
             str_replace_all("'", "")) #change name
  return(d)
}

#run function
#out <- lapply(studies, get_hits)
dm <- do.call(rbind, out)

#save all SNPs 
all_snps <- dm %>%
 select(rsid) %>%
 distinct()
write.table(all_snps, "/g/steinmetz/project/otar/data/gwas_selection/snplist_all_GWAS_traits.tsv", quote = F, sep ="\t", col.names = F, row.names = F)

#are all GWAS SNPs we included in here?
table(snps$rsid %in% all_snps$rsid) #no because some SNPs were included through CHEERS and they were not tested for the GWAS

#Annoyingly, some of these SNPs have two positions in this dataframe. From manual checks it seems *usually* the upstream variant is the correct one
#For plotting it does not make a huge difference, but it is important that we use the same position for the same SNP if it's found in multiple diseases
dm <- dm %>% 
  group_by(rsid) %>%
  mutate(position_updated = min(position)) %>% #keep just one (the lowest) position per rs ID
  ungroup()

#Perform liftover for all variants to get hg38 locations
# chain file for hg19 to hg38
#download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz', 'hg19ToHg38.over.chain.gz')
#R.utils::gunzip('hg19ToHg38.over.chain.gz')
chainObject <- import.chain('hg19ToHg38.over.chain')

#make GRanges object out of the GWAS SNPs
grObject <- GRanges(
  seqnames = paste0("chr", dm$chr),
  ranges = IRanges(start = dm$position_updated, end = dm$position_updated),
  rsid = dm$rsid)

#perform liftover
all_snps_hg38 <- as.data.frame(liftOver(grObject, chainObject))[,c('seqnames','start','end', 'rsid')] %>%
  distinct()

#Combine with dataframe to add new position
dm_hg38 <- dm %>%
  inner_join(all_snps_hg38, by = "rsid") %>%
  select(rsid, chr.hg38 = seqnames, pos.hg38 = start, p, beta, se, trait) %>%
  mutate(selected = case_when(rsid %in% snps$rsid ~ trait,
                              .default = "Z"))

#Save this information
write.table(dm_hg38, "/g/steinmetz/project/otar/data/gwas_selection/snplist_all_GWAS_traits_info_poshg38.tsv", quote = F, sep ="\t", col.names = T, row.names = F)

#### START HERE FOR PLOTTING ####
dm_hg38 <- fread(paste0(root, "data/gwas_selection/snplist_all_GWAS_traits_info_poshg38.tsv"))

#filter out HLA region chr6:28,510,120-33,480,577 for the plot
dm_hg38 <- dm_hg38 %>%
  filter(!(chr.hg38 == "chr6" & pos.hg38 > 28510120 & pos.hg38 < 33480577)) 

#for the promoters, remove AL and AC genes (which we also did for ordering)
proms_list <- proms %>% 
  filter(!str_detect(Gene, "^(AL|AC)[0-9]+\\.[0-9]+")) %>% #remove the AL and AC genes that remained
  select(gene_id)

#annotate gene expression table with genes that were included for one reason or another
genes_exp <- genes_exp %>%
  mutate(selection = case_when(gene_id %in% genes$gene_id ~ "enhancer target",
                               gene_id %in% tf$ENSEMBL ~ "transcription factor",
                               gene_id %in% proms_list$gene_id ~ "regular",
                               .default = "not selected"))

#for the enhancers, add activity
enh_act <- enh_act %>%
  mutate(peak = paste0(chr, ":", start+1000, "-", end-1000), #make the peaks without the extended windows to match enh list
         start.adj = start+1000,
         end.adj = end-1000) %>% 
  select(chr, start.adj, end.adj, peak, starts_with("ATAC"))

enh_act_comb <- enhs %>%
  full_join(enh_act, by = c("chr", "peak", "start" = "start.adj", "end" = "end.adj")) %>% 
  mutate(selected = case_when(str_detect(trait, "control") ~ "control",
                              str_detect(trait, "validation") ~ "control",
                              str_detect(trait, "eQTL") ~ "eQTL",
                              str_detect(trait, "GRN") ~ "GRN",
                              !is.na(trait) ~ "disease",
                              is.na(trait) ~ "Z"))
#                              .default = trait))


#### PLOT ####

#GWAS MANHATTAN

#ordering chromosomes
chrs <- paste0("chr", c(as.character(1:22), "X", "Y"))

#colour palette
palette <- c("Z" = "grey90",
             "Allergic_disease" =  "#A1D99B", #greens
             "Asthma" = "#41AB5D",
             "Eczema" = "#238B45",
             "Psoriasis" = "#006D2C",
             
             "Celiac_disease" = "#FDD0A2", #oranges
             "Crohns_disease" = "#FDAE6B",
             "Inflammatory_bowel_disease" = "#FD8D3C",
             "Ulcerative_colitis" = "#D94801",

             "Ankylosing_spondylitis" = "#9ECAE1", #blues
             "Gout" = "#6BAED6",
             "Rheumatoid_arthritis" = "#4292C6",
             "Systemic_lupus_erythematosus" = "#2171B5",
             "Multiple_sclerosis" = "#084594",
             
             "Primary_biliary_cholangitis" = "#9E9AC8", #purples
             "Primary_biliary_cirrhosis" = "#807DBA",
             "Primary_sclerosing_cholangitis" = "#6A51A3",
             
             "Type_1_diabetes" = "#FED976") #yellow 

# Compute (cumulative) chromosome sizes
chr_lengths <- dm_hg38 %>% 
  mutate(chr.hg38 = fct_relevel(chr.hg38, chrs)) %>%
  arrange(chr.hg38, pos.hg38) %>%
  group_by(chr.hg38) %>% 
  summarise(chr_len=as.numeric(max(pos.hg38, na.rm = T))) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len)

#Make plot data
snp_plot <- dm_hg38 %>% 
  
  # Add this info to the initial dataset
  left_join(chr_lengths, by="chr.hg38") %>%
  
  # Add a cumulative position of each SNP
  mutate(BPcum=pos.hg38+tot) 

#Make axis data frame
axisdf = snp_plot %>%
  group_by(chr.hg38) %>%
  mutate(chr.hg38 = str_replace(chr.hg38, "chr", "")) %>%
  summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

#Plot
p <- ggplot(snp_plot, aes(x=BPcum, y=-log10(p))) +
  
  #Plot all background points
  geom_point(data = snp_plot[snp_plot$selected == "Z",], 
             colour = "grey90", 
             alpha=0.1, 
             size = 1) +
  #Plot all foreground points
  geom_point(data = snp_plot[snp_plot$selected != "Z",], 
             aes(color=as.factor(selected)), #colour by trait
             alpha=0.8, #alpha different for highlighted vs. background
             size = 1) +
  scale_color_manual(values = palette) +
  #scale_alpha_manual(values = c(rep(0.8, 17), 0.1)) +
  
  # custom X axis:
  scale_x_continuous(limits = c(min(chr_lengths$tot), max(chr_lengths$tot)), label = axisdf$chr.hg38, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  
  #Change x and y text
  xlab("\nChromosome") +
  ylab("-log10(p)\n") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#rasterize the plot to make it much smaller
p.rast <- rasterize(p, dpi = 300)

#save this plot
pdf(paste0(outdir, "Manhattan_GWAS.pdf"), height = 4, width = 16)
p.rast
dev.off()

#ENHANCER MANHATTAN

#colour palette
enh_palette <- c("GRN" = "#D95F02", #orange
                 "eQTL" = "#1B9E77", #green
                 "disease" = "#7570B3", #purple
                 "control" = "grey90",
                 "Z" = "grey90") #background

#Make plot data
enh_plot <- enh_act_comb %>% 
  mutate(chr.hg38 = fct_relevel(chr, chrs), pos.hg38 = start) %>%
  arrange(chr.hg38, pos.hg38) %>%
  left_join(chr_lengths, by="chr.hg38") %>%
  # Add a cumulative position of each SNP
  mutate( BPcum=pos.hg38+tot) 

#Plot
#p_e <- ggplot(enh_plot, aes(x=BPcum, y=(ATAC_naive_16H_TH0.txt-ATAC_naive_16H_UNS.txt))) +
p_e <- ggplot(enh_plot, aes(x=BPcum, y=ATAC_naive_16H_TH0.txt)) +
  
#  geom_point(aes(color=as.factor(selected), #colour by trait
#                alpha=selected), #alpha different for highlighted vs. background
#            size = 1) +
  geom_point(data = enh_plot[enh_plot$selected == "Z",], 
             colour = "grey90", 
             alpha=0.1, 
             size = 1) +
  geom_point(data = enh_plot[enh_plot$selected != "Z",], 
             aes(color=as.factor(selected)), #colour by trait
             alpha=0.8, #alpha different for highlighted vs. background
             size = 1) +
  scale_color_manual(values = enh_palette) +
#  scale_alpha_manual(values = c(rep(0.8, 4), 0.05)) +
  
  # custom X axis:
  scale_x_continuous(limits = c(min(chr_lengths$tot), max(chr_lengths$tot)), label = axisdf$chr.hg38, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +     # remove space between plot area and x axis
  
  #Change x and y text
  xlab("\nChromosome") +
  ylab("Activity in naive Th0 T-cells\n") +
#  ylab("Activity in Th0 - activity in unstimulated naive T-cells\n") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#rasterize the plot to make it much smaller
p_e.rast <- rasterize(p_e, dpi = 300)

#save this plot
pdf(paste0(outdir, "Manhattan_enhancers.pdf"), height = 4, width = 16)
p_e.rast
dev.off()

#PROMOTER MANHATTAN

#colour palette
prom_palette <- c("regular" = "#386CB0",
             "enhancer target" =  "#7FC97F",
             "transcription factor" = "#F0027F")

#Make plot data
prom_plot <- genes_exp %>% 
  mutate(chr.hg38 = fct_relevel(chr, chrs), pos.hg38 = start) %>%
  arrange(chr.hg38, pos.hg38) %>%
  left_join(chr_lengths, by="chr.hg38") %>%
  # Add a cumulative position of each SNP
  mutate( BPcum=pos.hg38+tot) 

#Plot
p_p <- ggplot(prom_plot, aes(x=BPcum, y=tpm)) +
#  geom_point(aes(color=selection), #colour by whether they are an enhancer screen target, TF or regular promoter
#              size = 1, alpha = 0.2) +
  geom_point(data = prom_plot[prom_plot$selection == "not selected",], 
             colour = "grey90",
             size = 1, 
             alpha = 0.1) +
  geom_point(data = prom_plot[prom_plot$selection != "not selected",],
             aes(color=selection), #colour by whether they are an enhancer screen target, TF or regular promoter
             size = 1, 
             alpha = 0.8) +
  scale_color_manual(values = prom_palette) +
  scale_x_continuous(limits = c(min(chr_lengths$tot), max(chr_lengths$tot)), label = axisdf$chr.hg38, breaks= axisdf$center) +
  scale_y_log10() +
  xlab("\nChromosome") +
  ylab("Gene expression level (TPM)\n") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#rasterize the plot to make it much smaller
p_p.rast <- rasterize(p_p, dpi = 300)

#save this plot
pdf(paste0(outdir, "Manhattan_promoters.pdf"), height = 4, width = 16)
p_p.rast
dev.off()


#COMBINE MANHATTAN PLOTS
#combined <- p.rast / p_e.rast / p_p.rast
combined <- p.rast / p_e.rast / p_p.rast + plot_layout(guides = "collect") & theme(axis.title.x=element_blank()) #to ensure they use the same axis

#save this plot
pdf(paste0(outdir, "Manhattan_combined.pdf"), height = 8, width = 16)
combined
dev.off()


#SNP STACKED BARPLOT

#order of diseases for plot
levels_diseases <- c("Allergic_disease", "Asthma", "Eczema", "Psoriasis", "Skin and allergic",
                     "Celiac_disease", "Crohns_disease", "Inflammatory_bowel_disease", "Ulcerative_colitis", "Gastrointestinal",
                     "Ankylosing_spondylitis", "Gout", "Rheumatoid_arthritis", "Systemic_lupus_erythematosus", "Multiple_sclerosis", "Rheumatic",
                     "Primary_biliary_cholangitis", "Primary_biliary_cirrhosis", "Primary_sclerosing_cholangitis", "Hepatic",
                     "Type_1_diabetes", "Multiple diseases")

#get data formatted for the plot
barplot_data <- snps %>%
  filter(id != "nonGWAS") %>%
  select(rsid, trait, disease_group) %>%
  distinct() %>%
  group_by(rsid) %>%
  mutate(disease_count = n_distinct(trait),                 # Count unique diseases linked to the SNP
         disease_group_count = n_distinct(disease_group),   # Count unique disease groups linked to the SNP
         condition = case_when(
           disease_count == 1 ~ trait,                     # Case 1: SNP linked to only one disease
           disease_group_count == 1 ~ disease_group,       # Case 2: SNP linked to multiple diseases from one group
           TRUE ~ "Multiple diseases"                      # Case 3: SNP linked to multiple diseases from different groups
         )) %>%
  ungroup() %>% # Ungroup when done
  mutate(#condition = str_replace_all(condition, "_", " "),
         condition = fct_relevel(condition, levels_diseases),
         dummy = 1) %>%
  select(rsid, condition, dummy) %>%
  distinct()

#add overarching groups to colour palette
palette_bar <- c(palette, 
                 "Skin and allergic" = "#00441B",
                 "Gastrointestinal" = "#A63603",
                 "Rheumatic" = "#08306B",
                 "Hepatic" = "#3F007D",
                 "Multiple diseases" = "grey20")

# Summarize data to get the counts and positions for labels
label_data <- barplot_data %>%
  group_by(condition) %>%
  summarise(count = dplyr::n()) %>%  # Explicitly use dplyr::n() to count rows
  arrange(match(condition, rev(levels_diseases))) %>%
  mutate(cumsum = cumsum(count) - (0.5 * count))  # Calculate center of each stacked bar

#bar plot
p <- ggplot(barplot_data) +
  geom_bar(aes(x = dummy, fill = condition), stat = "count", position = "stack") +
  scale_fill_manual(values = palette_bar) +
  geom_text(data = label_data, aes(x = 1.05, y = cumsum, label = condition), 
            hjust = 0, size = 5) +  # Adjust hjust to left-align
  coord_cartesian(xlim = c(0.8, 1.2)) +  # Adjust x-axis limits to make space for text
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

#save this plot
pdf(paste0(outdir, "Barplot_GWAS.pdf"), height = 10, width = 8)
p
dev.off()


 #OLD CODE 


#Get the GWAS results (with p-values) for the SNPs we ended up including
# get_hits <- function(study){
#   cat("Now submitting", study, "\n")
#   d <- tophits(id=study,
#                 opengwas_jwt = opengwas_jwt)
#   cat(length(unique(d$rsid)), "SNPs found\n")
#   d <- d %>%
#     mutate(trait = trait %>%
#              str_replace_all(" \\(.*", "") %>%
#              str_replace_all(" ", "_") %>%
#              str_replace_all("'", ""))
#   return(d)
# }
# 
# retrieve_snp_info <- function(study){
#   disease <- snps %>%
#     filter(id == !!study) %>%
#     select(trait) %>%
#     distinct() %>%
#     pull(trait)
#   cat("\nNow submitting", study, "=", disease, "\n")
#   snps_disease <- snps %>%
#     filter(id == !!study) %>%
#     pull(rsid)
#   cat(length(unique(snps_disease)), "SNPs in query list\n")
#   d <- associations(variants = snps_disease, 
#                     id = study,
#                     opengwas_jwt = opengwas_jwt)
#   cat(length(unique(d$rsid)), "SNPs found in database\n")
#   return(d)
# }
# 
# #run function
# #out <- lapply(studies, retrieve_snp_info)
# dm <- do.call(rbind, out)
# 
# length(unique(dm$rsid))
# length(unique(snps$rsid))
# #issue: not all SNPs in our original list have a location here
# 
# #combine output with original SNP list to get their Hg38 position
# snp_info <- snps %>%
#   left_join(dm, by = c("rsid", "id")) %>%
#   select(rsid, chr.hg38 = chr.x, pos.hg38 = start,
#          chr.hg19 = chr.y, pos.hg19 = position,
#          trait = trait.x, study.id = id,
#          se, p, beta, n, ea, nea)
# length(unique(snp_info$rsid))
# #issue: not all SNPs in our original list have effect sizes here
# 
# #probably due to proxies that I got back then for each SNP



#THIS PART GAVE EMPTY LINES IN LIFTOVER!!!
# #save all SNPs per chromosome so it's easier to get the positions
# chrs <- c(as.character(1:22), "X", "Y")
# 
# walk(chrs, function(chr) {
#   chr_snps <- dm %>%
#     filter(chr == !!chr) %>%
#     select(rsid) %>%
#     distinct()
#   
#   # Define the output file name
#   file_name <- paste0("/g/steinmetz/project/otar/data/gwas_selection/snplist_all_GWAS_traits_chr", chr, ".tsv")
#   
#   # Write the SNPs for this chromosome to a file
#   write.table(chr_snps, file_name, quote = F, sep = "\t", col.names = F, row.names = F)
# })
# #Perform liftover for all variants to get hg38 locations, split by chr
# #script: /g/steinmetz/project/otar/scripts/get_hg38_coordinates_for_all_gwas_snps.sh
# all_snps_hg38 <- fread("/g/steinmetz/project/otar/data/gwas_selection/snplist_all_GWAS_traits_poshg38.tsv")
# files_hg38 <- list.files(path = "/g/steinmetz/project/otar/data/gwas_selection/", pattern = "_poshg38.tsv", full.names = T)
# 
# all_snps_hg38 <- do.call(rbind, lapply(files_hg38, fread))
