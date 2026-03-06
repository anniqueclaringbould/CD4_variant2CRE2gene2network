#############################################
## Script to quantify enrichment of PCHi-C ##
########### DpnII-fragment level ############
#############################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(plyranges)

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
hicdir <- paste0(root, "interpretation/hi-c_enrichment/")
indir <- paste0(hicdir, "data/")
outdir <- paste0(hicdir, "output/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "design/gene_selection/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read screen results, combined and annotated
res <- fread(paste0(root, "interpretation/chromatin_analyses/results/results_df_with_promoterC_annotated.csv"))

#read matrix with significant PCHi-C results
pchic <- fread(paste0(indir, "CD4_chicago_fres_5kb_abc_023_fres_extended_peakm_13012025.txt"))

#read all tested promoters (baits) from PCHi-C experiment
promoters_dpn <- fread(paste0(indir, "hg38_dpnII.baitmap.headers.txt"))
promoters_5kb <- fread(paste0(indir, "human_DpnII_5K_sol_baits.baitmap.headers.txt"))

#read all tested DnpII fragments from PCHi-C experiment
regions_dpn <- fread(paste0(indir, "hg38_dpnII.rmap.headers.txt"))
regions_5kb <- fread(paste0(indir, "human_DpnII_5K_sol_baits.rmap.headers.txt"))

#### ADJUST ####

#Plot fragment lengths in Dpn-II resolution
regions_dpn %>% mutate(fragment_length = oe_end_fres - oe_start_fres) %>% 
  ggplot(., aes(x = fragment_length)) + 
  geom_histogram(bins = 100) + 
  xlim(0,2000) + 
  theme_classic() +
  ggtitle("Distribution of DpnII fragment lengths")

regions_dpn %>% mutate(fragment_length = oe_end_fres - oe_start_fres) %>% summary()

# for control perturbations which were in multiple panels, select the most significant one
res.ctrl <- res %>% 
  filter(ctrl_perturbation == TRUE) %>% 
  group_by(response_id, grna_target) %>% 
  slice_min(p_value, n = 1, with_ties = FALSE)

# add back to results table with discovery pairs
res <- res %>% 
  filter(ctrl_perturbation == FALSE) %>% 
  bind_rows(res.ctrl) %>% 
  arrange(desc(p_value)) 

#Filter the results to test
res.filt <- res %>%
  filter(pert_chr == gene_chr) %>% #keep enhancer-gene links on the same chromosome
  filter(pass_qc == TRUE) %>% #keep enhancer-gene links that were tested and pass QC
  filter(class_manuscript == "intergenic-distal" | true_enhancer == TRUE) %>%
  mutate(link_id = paste0(grna_target, "_", gene_id)) 

#Adjust PCHi-C results at fragment resolution
pchic_dpn <- pchic %>%
  mutate(gene = baitName) %>%
  separate_rows(gene, sep = "/") %>% #some bait names have multiple genes, split to one row per gene
  filter(!is.na(chicago_score_fres)) %>% #keep only fragment-level results
  select(baitChr, bait_start_fres, bait_end_fres, baitID_fres, baitName, gene, #keep only fragment-related columns
         oeChr, oe_start_fres, oe_end_fres, oeID_fres,
         dist, N_fres, chicago_score_fres) %>%
  distinct()
#no missing bait names or anything else

#Adjust PCHi-C results at 5kb resolution
#Some entries have empty values for oeChr (element chromosome), even if their region IDs are present
#Some entries have empty baitNames (promoter names), even if their promoter IDs do have one associated
#Therefore, left_join the regions and promoters with their information based on element ID
pchic_5kb <- pchic %>%
  select(-oeChr) %>%
  left_join(regions_5kb, by = c("oeID_5kb")) %>%
  left_join(promoters_5kb, by = c("baitID_5kb")) %>%
  mutate(baitName = case_when(baitName.x != "" ~ baitName.x, #if it's not missing, use original baitname
                              baitName.x == "" & !is.na(baitName.y) ~ baitName.y, #if it's missing but there is a value in the annotation, use annotation
                              baitName.x == "" & is.na(baitName.y) ~ baitName.y), #if it's missing and there is nothing in the annotation, keep in as NA value
         gene = baitName) %>%
  separate_rows(gene, sep = "/") %>% #some bait names have multiple genes, split to one row per gene
  filter(!(is.na(chicago_score_5kb) & is.na(ABC.Score))) %>% #keep only fragment-level results
  select(baitChr = baitChr.x, bait_start_5kb, bait_end_5kb, baitID_5kb, baitName, gene, #keep only 5kb-related columns
         oeChr, oe_start_5kb, oe_end_5kb, oeID_5kb,
         dist, N_5kb, N_abc, chicago_score_5kb, ABC.Score) %>%
  distinct()

### Overlap tested bait (promoter) regions from PCHi-C data with targets tested in our data ###

#Select genes from our results
genes <- res.filt %>%
  select(starts_with("gene"), response_id) %>%
  distinct()

#Select genes from promoter baits
promoter_genes_dpn <- promoters_dpn %>%
  select(baitName) %>%
  mutate(gene = baitName) %>%
  separate_rows(gene, sep = "/") %>% #some bait names have multiple genes here
  distinct()

promoter_genes_5kb <- promoters_5kb %>%
  select(baitName) %>%
  mutate(gene = baitName) %>%
  separate_rows(gene, sep = "/") %>% #some bait names have multiple genes here
  distinct()

#Number of overlapping genes
length(intersect(genes$response_id, promoter_genes_dpn$gene))
length(intersect(genes$response_id, promoter_genes_5kb$gene)) #same number of genes, same number of elements too
table(promoter_genes_dpn$gene == promoter_genes_5kb$gene) #exactly the same genes, so can use 1 list for the gene overlap

#Filter our results for list of genes that were also tested in PCHi-C experiment
res_filtered_genes <- res.filt %>%
  filter(response_id %in% promoter_genes_dpn$gene)

#Number of genes in results dataframe
length(unique(res_filtered_genes$response_id)) #1413 genes


### Overlap tested fragments with enhancers tested in our data ###

#Select regulatory elements from our results
cres <- res_filtered_genes %>%
  select(grna_target, pert_chr, pert_start, pert_end, ctrl_perturbation, element_type) %>%
  distinct()

#make CREs into GRanges object 
gr.cres <- cres %>% 
  makeGRangesFromDataFrame(., seqnames.field="pert_chr", start.field="pert_start", end.field="pert_end", keep.extra.columns=T)
seqlevelsStyle(gr.cres) <- "UCSC" #consistent gene formatting

#make PCHi-C dpnII regions info into GRanges object
gr.regions_dpn <- regions_dpn %>%
  makeGRangesFromDataFrame(., seqnames.field="oeChr", start.field="oe_start_fres", end.field="oe_end_fres", keep.extra.columns=T)
seqlevelsStyle(gr.regions_dpn) <- "UCSC" #consistent gene formatting

#make PCHi-C 5kb regions info into GRanges object
gr.regions_5kb <- regions_5kb %>%
  makeGRangesFromDataFrame(., seqnames.field="oeChr", start.field="oe_start_5kb", end.field="oe_end_5kb", keep.extra.columns=T)
seqlevelsStyle(gr.regions_5kb) <- "UCSC" #consistent gene formatting

#overlap dpnII regions based on location 
gr.regions_overlap_dpn <- join_overlap_inner(gr.cres, gr.regions_dpn, maxgap = 0)
length(unique(gr.regions_overlap_dpn$grna_target)) == length(unique(res_filtered_genes$grna_target)) 
#307 regions

#overlap 5kb regions based on location 
gr.regions_overlap_5kb <- join_overlap_inner(gr.cres, gr.regions_5kb, maxgap = 0)
length(unique(gr.regions_overlap_5kb$grna_target)) == length(unique(res_filtered_genes$grna_target)) 
#307 regions

#No need to further filter our results for list of regions that were also tested in PCHi-C experiment
#because all our CREs have a match in PCHi-C experiment


### Join our results and PCHi-C results by merging CREs/regions then filtering for same target gene ###

#make results into GRanges object based on the location of the CRE 
gr.res <- res_filtered_genes %>% 
  makeGRangesFromDataFrame(., seqnames.field="pert_chr", start.field="pert_start", end.field="pert_end", keep.extra.columns=T)
seqlevelsStyle(gr.res) <- "UCSC" #consistent gene formatting
  
#make PCHi-C dpnII results into GRanges object based on the location of the region (not the promoter gene)
gr.pchic_dpn <- pchic_dpn %>%
  makeGRangesFromDataFrame(., seqnames.field="oeChr", start.field="oe_start_fres", end.field="oe_end_fres", keep.extra.columns=T)
seqlevelsStyle(gr.pchic_dpn) <- "UCSC" #consistent gene formatting

#make PCHi-C 5kb results into GRanges object based on the location of the region (not the promoter gene)
gr.pchic_5kb <- pchic_5kb %>%
  makeGRangesFromDataFrame(., seqnames.field="oeChr", start.field="oe_start_5kb", end.field="oe_end_5kb", keep.extra.columns=T)
seqlevelsStyle(gr.pchic_5kb) <- "UCSC" #consistent gene formatting

#overlap significant PCHi-C links with our results (all tested, not just significant), based on region
gr.res_overlap_dpn <- join_overlap_inner(gr.res, gr.pchic_dpn, maxgap = 0) 
length(unique(gr.res_overlap_dpn$grna_target)) #only 35 overlap

gr.res_overlap_5kb <- join_overlap_inner(gr.res, gr.pchic_5kb, maxgap = 0) 
length(unique(gr.res_overlap_5kb$grna_target)) #only 259 overlap

#the regions-only overlap (1036 regions) was on all tested regions from the rmap
#this overlap is based on the PCHi-C results file, which only includes significant results

#filter to keep only results where the PCHi-C promoter gene also matches the target gene in our results
res_overlap_filtered_dpn <- gr.res_overlap_dpn %>%
  as.data.frame() %>%
  filter(response_id == gene)

res_overlap_filtered_5kb <- gr.res_overlap_5kb %>%
  as.data.frame() %>%
  filter(response_id == gene)

#make a list of results that overlap an enhancer screen effect
examples_dpn <- res_overlap_filtered_dpn %>%
  filter(baitName != "",
         significant == TRUE) %>%
  select(link_id, response_id, is_target, ctrl_perturbation, element_type, class_manuscript, 
         baitID_fres, baitName, oeID_fres, chicago_score_fres) %>%
  distinct()
fwrite(examples_dpn, paste0(outdir, "ETPs_significant_validated_by_PCHi-C_DnpII_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

examples_5kb <- res_overlap_filtered_5kb %>%
  filter(baitName != "",
         significant == TRUE) %>%
  select(link_id, response_id, is_target, ctrl_perturbation, element_type, class_manuscript, 
         baitID_5kb, baitName, oeID_5kb, chicago_score_5kb, ABC.Score) %>%
  distinct()
fwrite(examples_5kb, paste0(outdir, "ETPs_significant_validated_by_PCHi-C_5kb_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

#types of interactions being substantiated by these links
examples_dpn %>%
  select(link_id, class_manuscript) %>%
  distinct() %>%
  group_by(class_manuscript) %>%
  count()

examples_5kb %>%
  select(link_id, class_manuscript) %>%
  distinct() %>%
  group_by(class_manuscript) %>%
  count()
  
#select only link ID, significance, and bait info for dpnII resolution
res_overlap_filtered_dpn2 <- res_overlap_filtered_dpn %>%
  select(link_id, significant, baitName, baitID_fres) %>%
  dplyr::rename(baitID = baitID_fres) %>%
  mutate(resolution = "dpnII")

#select only link ID, significance, and bait info for 5kb resolution
res_overlap_filtered_5kb2 <- res_overlap_filtered_5kb %>%
  select(link_id, significant, baitName, baitID_5kb) %>%
  dplyr::rename(baitID = baitID_5kb) %>%
  mutate(resolution = "5kb")

#consolidate overlapping results for both resolutions to one line per peak-gene combination
res_overlap_minimal <- bind_rows(res_overlap_filtered_dpn2, res_overlap_filtered_5kb2) %>%
  distinct() %>%
  group_by(link_id, significant, resolution) %>%
  summarise(
    promoter_names   = paste(unique(baitName), collapse = ","),
    promoter_baitIDs = paste(unique(baitID), collapse = ","),
    .groups = "drop"
  ) %>%
  group_by(link_id, significant) %>%
  summarise(
    pchic_results = paste0("res ", resolution, ": ", promoter_baitIDs, collapse = "; "),
    promoter_names = paste(unique(promoter_names), collapse = "; "),
    .groups = "drop"
  )

#annotate general results with PCHi-C overlaps, keep one line per enh-gene link
res_pchic_ann_minimal <- res_filtered_genes %>%
  left_join(res_overlap_minimal, by = c("link_id", "significant")) %>%
  select(link_id, response_id, is_target, class_manuscript, pchic_results, promoter_names, significant) %>%
  mutate(pchic_overlap = ifelse(!is.na(promoter_names), TRUE, FALSE))
fwrite(res_pchic_ann_minimal, paste0(outdir, "Screen_results_minimal_annotated_PCHi-C_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

#total number of genes
res_pchic_ann_minimal %>% filter(significant == T) %>% ungroup() %>% filter(pchic_overlap == T) %>% select(response_id) %>% distinct() %>% nrow()

#significant vs. non-significant effects
cont.table.sign <- table(res_pchic_ann_minimal$pchic_overlap,
                         res_pchic_ann_minimal$significant)

#Fisher test
fisher.test(cont.table.sign)
pval <- fisher.test(cont.table.sign)$p.value
odds_ratio <- fisher.test(cont.table.sign)$estimate
log_odds_ratio <- log(odds_ratio)

#Calculate proportions of significant results
proportion_data <- res_pchic_ann_minimal %>%
  group_by(pchic_overlap) %>%
  filter(!is.na(pchic_overlap)) %>%
  summarize(
    count_significant = sum(significant == TRUE),     # Count of significant results
    total_count = dplyr::n(),                                  # Total count of results
    pct_significant = count_significant / total_count * 100 # Calculate proportion
  ) %>%
  mutate(pchic = case_when(pchic_overlap == T ~ "pcHi-C overlap",
                            pchic_overlap == F ~ "no pcHi-C overlap"))

# Bar plot of proportions
y_max <- max(proportion_data$pct_significant)

ggplot(proportion_data, aes(x = pchic, y = pct_significant, fill = pchic)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("grey70", "#0263a3")) +
  geom_segment(
    aes(x = 1, xend = 2, y = y_max * 1.05, yend = y_max * 1.05),
    inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = y_max * 1.12, label = paste0("OR = ", round(odds_ratio, digits = 1), "\np = ", signif(pval, 3))) +
  labs(y = "Percentage of results that is significant (%)\n") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) 

ggsave(paste0(outdir, "Pct_significant_results_pcHi-C_", Sys.Date(), ".pdf"), height = 5, width = 3.5)
  
res_pchic_ann_minimal %>%
  mutate(pchic = case_when(pchic_overlap == T ~ "pcHi-C\noverlap",
                            pchic_overlap == F ~ "no pcHi-C\noverlap")) %>% 
  ggplot(.) +
  geom_bar(aes(x = pchic, fill = significant), stat = "count", position = "fill", colour = "black") + 
  scale_fill_manual(values = c("grey70", "#0263a3")) +
  geom_segment(
    aes(x = 1, xend = 2, y = 1.05, yend = 1.05),
    inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = 1.12, label = paste0("OR = ", round(odds_ratio, digits = 1), "\np = ", signif(pval, 3))) +
  labs(y = "Proportion of results\n", x = "") +
  theme_bw()

ggsave(paste0(outdir, "Prop_significant_results_pcHi-C_", Sys.Date(), ".pdf"), height = 5, width = 3.5)










# #OLD
# #make genes into GRanges object 
# gr.genes <- genes %>% 
#   makeGRangesFromDataFrame(., seqnames.field="gene_chr", start.field="gene_start", end.field="gene_end", keep.extra.columns=T)
# seqlevelsStyle(gr.genes) <- "UCSC" #consistent gene formatting
# 
# #make promoter info into GRanges object
# gr.promoters <- promoters %>%
#   makeGRangesFromDataFrame(., seqnames.field="baitChr", start.field="bait_start_fres", end.field="bait_end_fres", keep.extra.columns=T)
# seqlevelsStyle(gr.promoters) <- "UCSC" #consistent gene formatting
# 
# #examples of doubles
# tmp <- gr.genes_overlap %>%
#   as.data.frame() %>%
#   select(gene_id, baitName, response_id) %>%
#   distinct() 
# 
# double_genes <- tmp[duplicated(tmp$gene_id),]
# 
# genes %>% filter(gene_chr == "chr11") %>% filter(gene_id == "ENSG00000002330")
# promoters %>% filter(baitChr == 11) %>% filter(bait_start_fres > 64275483) %>% filter(bait_end_fres < 64287078)
# 
# genes %>% filter(gene_id == "ENSG00000010610")
# promoters %>% filter(baitChr == 12) %>% filter(bait_start_fres > 6780000) %>% filter(bait_start_fres < 6821000)
# 
# gr.genes_overlap %>%
#   as.data.frame() %>%
#   filter(gene_id == "ENSG00000010610")



#####LEFT HERE: CHECK IF ABOVE IS CORRECT
#table(pchic_5kb$baitChr.x == pchic_5kb$baitChr.y, useNA = "always")
#table(pchic_5kb$baitName.x == pchic_5kb$baitName.y, useNA = "always")



# pchic3 <- pchic2 %>%
#   left_join(promoters, by = c("baitID_fres", "bait_start_fres", "bait_end_fres"))
# 
# #!!!!!! ISSUES
# pchic3 %>% 
#   filter(baitName.x == "") %>% 
#   filter(!is.na(baitName.y)) %>%
#   select(baitID_fres) %>%
#   distinct() %>%
#   nrow()
# #~240k results have no baitname in PCHi-C results, but do have a baitName in rmap (example: ID 1970-NOC2L, 1978-PLEKHN1)
# #8508 baits
# 
# pchic3 %>% 
#   filter(baitName.x == "") %>% 
#   filter(is.na(baitName.y)) %>%
#   select(baitID_fres) %>%
#   distinct() %>%
#   nrow()
# #~790k results have no baitname in PCHi-C results, and also not in rmap (example: ID 2, 3, 4, 5..)
# #29014 baits
# 
# pchic3 <- pchic2 %>%
#   select(-baitName) %>%
#   left_join(promoters, by = c("baitID_fres", "bait_start_fres", "bait_end_fres"))
# 
