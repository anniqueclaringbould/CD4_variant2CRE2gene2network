########################################
#### Script to plot eQTL enrichment ####
## eQTL catalog + eQTLGen ##############
########################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
outdir <- paste0(root, "interpretation/eqtl_enrichment/output/")

#### READ DATA ####

#read eQTL enrichment from eQTLGen p2
res_eQTL_eqtlgen_ann_minimal <- fread(paste0(outdir, "Screen_results_minimal_annotated_eQTL_2025-12-23.tsv"))

#read eQTL enrichment from eQTL catalog datasets

#list of files
files <- list.files(path = paste0(outdir, "eQTL_catalogue/"), pattern = "Screen")

#function to read
read_one <- function(f, indir) {
  dt <- fread(file.path(indir, f))
  
  dataset <- sub("^Screen_results_minimal_annotated_eQTL_", "", f)
  dataset <- sub("\\.tsv$", "", dataset)
  
  dt_out <- dt %>%
    mutate(!!paste0("eQTL_", dataset) := !is.na(variant) & variant != "") %>%
    select(link_id, response_id, is_target, interaction_type, significant, !!paste0("eQTL_", dataset))
  
  #setnames(dt_out, "overlaps_eQTL", paste0("eQTL_", dataset))
  dt_out
}

#apply function
lst <- lapply(files, read_one, indir = paste0(outdir, "eQTL_catalogue/"))

#combine into one data frame
res_combined <- Reduce(function(x, y)
  merge(
    x, y,
    by = c("link_id", "response_id", "is_target", "interaction_type", "significant"),
    all = TRUE
  ),
  lst
)

#make one column with yes/no eQTL and another containing the dataset(s)
eqtl_cols <- grep("^eQTL_", names(res_combined), value = TRUE)

res_combined2 <- res_combined %>%
  rowwise() %>%
  mutate(
    overlaps_eQTL_any = any(c_across(all_of(eqtl_cols))),
    overlaps_eQTL_datasets = if (overlaps_eQTL_any) {
      paste(eqtl_cols[c_across(all_of(eqtl_cols))], collapse = ";") %>%
        sub("^eQTL_", "", .)
    } else {
      NA_character_
    }
  ) %>%
  ungroup()

#Check for specific examples
res_combined2 %>%
  filter(str_detect(link_id, "chr16:10966558-10966877")) %>%
  filter(response_id == "DEXI") %>%
  select(overlaps_eQTL_any, overlaps_eQTL_datasets)



#INTERACTION TYPES
#compare interaction types between significant ETPs with eQTL and without eQTL
eqtl_interaction_types <- res_combined %>%
  filter(is_target == TRUE) %>% #keep only prioritised hits
  filter(significant == TRUE) %>% #keep only prioritised hits
  group_by(overlaps_eQTL_any, interaction_type) %>%
  tally() %>%
  mutate(eQTL = case_when(overlaps_eQTL_any == TRUE ~ "eQTL",
                          overlaps_eQTL_any == FALSE ~ "no eQTL")) %>%
  ungroup() %>%
  select(interaction_type, n, eQTL)

# Colour palette to match with prior colouring
# colors for each class
class_colors <- c(
  "promoter-proximal" = "#0263a3",   
  "promoter-distal" = "#4299d4", 
  "intragenic-proximal" ="#9ad8fc",
  "intragenic-distal" = "#e57373", 
  "intergenic-distal"= "#c62828" 
)

#Order of interaction types
interaction_levels <- c("promoter-proximal", "promoter-distal", "intragenic-proximal", "intragenic-distal", "intergenic-distal")

# Stacked bar plot of proportions
ggplot(eqtl_interaction_types, aes(fill=fct_rev(fct_relevel(interaction_type, interaction_levels)), y=n, x=fct_rev(eQTL))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = class_colors)  +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  ylab("\n% of significant prioritised enhancer-gene links") +
  xlab("") +
  guides(fill = guide_legend(title = "", reverse = T))  +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.text = element_text(size = 12))
ggsave(paste0(outdir, "Percentage_interaction_types_prioritised_effects_all_eQTLs_", Sys.Date(), ".pdf"), height = 2, width = 8)

# Stacked bar plot of absolute numbers
ggplot(eqtl_interaction_types, aes(fill=fct_rev(fct_relevel(interaction_type, interaction_levels)), y=n, x=fct_rev(eQTL))) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = class_colors)  +
  theme_classic() +
  coord_flip() +
  ylab("\nNumber of significant prioritised enhancer-gene links") +
  xlab("") +
  guides(fill = guide_legend(title = "", reverse = T))  +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.text = element_text(size = 12))
ggsave(paste0(outdir, "Number_interaction_types_prioritised_effects_all_eQTLs_", Sys.Date(), ".pdf"), height = 2, width = 8)




