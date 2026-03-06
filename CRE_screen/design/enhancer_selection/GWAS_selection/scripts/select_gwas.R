###########################################################
###### OTAR project: autoimmune trait GWAS selection ######
# Script to get SNP lists from a selection of GWAS traits #
###########################################################

#Input options
library(ieugwasr)
library(rsnps)
library(tidyverse)
library(data.table)

dir <- "/g/steinmetz/project/otar/data/gwas_selection"

#Get list of GWAS summary stats available on the OpenGWAS database
#https://gwas.mrcieu.ac.uk/
all_studies <- gwasinfo()

#Filter out the relevant studies
filtered_studies <- all_studies %>%
  filter(category == "Disease" | category == "NA") %>%
  filter(subcategory == "Autoimmune / inflammatory" | subcategory == "NA") %>%
  filter(population == "European" | population == "Mixed") %>% #select only studies done on European or mixed populations because that matches with our T cell data
  filter(sample_size > 10000) %>% #remove studies with fewer than 10,000 donors  
  filter(!str_detect(id, "eqtl")) %>% #remove any eQTL studies
  filter(str_detect(tolower(trait), "allergic disease") |
           tolower(trait) == "ankylosing spondylitis" |
           tolower(trait) == "asthma" |
           str_detect(tolower(trait), "celiac") |
           tolower(trait) == "crohn's disease" |
           tolower(trait) =="eczema" |
           tolower(trait) =="gout" |
           str_detect(tolower(trait), "inflammatory bowel") |
           str_detect(tolower(trait), "idiopathic arthritis") |
           str_detect(tolower(trait), "multiple sclerosis") |
           str_detect(tolower(trait), "primary biliary") |
           str_detect(tolower(trait), "primary sclerosing") |
           tolower(trait) == "psoriasis" |
           tolower(trait) == "rheumatoid arthritis" |
           str_detect(tolower(trait), "lupus") |
           str_detect(tolower(trait), "type 1 diabetes") |
           tolower(trait) == "ulcerative colitis") %>% #select autoimmune GWAS of interest
  mutate(trait = str_replace_all(trait, " \\(.*", "")) %>% #collapse trait versions by removing parentheses
  group_by(trait) %>%
  slice_max(ncase, with_ties = F) #keep study with highest number of cases

#Get GWAS hits per disease in the filtered list write out the rs IDs
#function to get hits, save & write out rs IDs
get_hits <- function(study){
  dm <- tophits(id=study)
  snps <- unique(dm$rsid)
  trait <- dm %>%
    select(trait) %>%
    distinct() %>%
    str_replace_all(., " \\(.*", "") %>%
    str_replace_all(., " ", "_") %>%
    str_replace_all(., "'", "")
  dir.create(paste0(dir, "/", trait))
  write.table(snps, paste0(dir, "/", trait, "/snplist.tsv"), quote = F, sep ="\t", col.names = F, row.names = F)
}

#make a list of studies and diseases
studies <- filtered_studies$id
diseases <- filtered_studies$trait %>%
  str_replace_all(., " ", "_") %>%
  str_replace_all(., "'", "")

#run function
lapply(studies, get_hits)

#save GWAS info and filtered study details
write.table(all_studies, paste0(dir, "/all_studies_OpenGWAS_", Sys.Date(), ".tsv"), quote = F, sep ="\t", col.names = T, row.names = F)
write.table(filtered_studies, paste0(dir, "/filtered_studies_OpenGWAS_", Sys.Date(), ".tsv"), quote = F, sep ="\t", col.names = T, row.names = F)
write.table(diseases, paste0(dir, "/disease_names.tsv"), quote = F, sep ="\t", col.names = F, row.names = F)
