###################################################################################################
#### Script to select enhancers that are more active in Th0 as compared to unstimulated states ####
######################## And extend the ATAC seq peaks by various windows #########################
###################################################################################################

#SET-UP
library(tidyverse)
library(data.table)

#dir <- "/g/steinmetz/project/otar/select_enhancer_peaks/normalised_peaks/"
dir <- "/g/steinmetz/project/otar/select_enhancer_peaks/normalised_peaks_across_conditions/"
source("/g/scb/zaugg/claringb/tools/personal_colour_palette.R") #get personal colour palette

peak <- "H3K27ac" #"ATAC" "H3K27ac"
peakname <- ifelse(peak == "ATAC", "ATAC", "ChM")
normalisation <- "normToMax_quantileNorm_euclideanNorm" #"normToMax" "normToMax_quantileNorm" "normToMax_quantileNorm_euclideanNorm"
file <- paste0(dir, peak, "_counts_", normalisation, ".txt") 

#READ DATA
dm <- fread(file)

#save column names
col_uns <- paste0(peakname, "_naive_16H_UNS.txt")
col_16h <- paste0(peakname, "_naive_16H_TH0.txt")
col_5d <- paste0(peakname, "_naive_D5_TH0.txt")

#PLOT
dm2 <- dm %>%
  select(chr, start, end, !!sym(col_16h), !!sym(col_5d), !!sym(col_uns)) %>%
  pivot_longer(cols = starts_with(peakname), names_to = "condition")

#plot density distribution of scores per condition of interest
if (normalisation == "normToMax_quantileNorm_euclideanNorm") {
  ggplot(dm2, aes(x = value)) +
    geom_density(aes(group = condition, col = condition, fill = condition), alpha = 0.1) +
    scale_color_manual(values = c(ga_midgreen, ga_navy, ga_darkgreen)) +
    scale_fill_manual(values = c(ga_midgreen, ga_navy, ga_darkgreen)) +
    theme_bw() +
    ggtitle(paste0("Peak activity \nNormalisation: ", normalisation))
} else {
  ggplot(dm2, aes(x = log10(value))) +
    geom_density(aes(group = condition, col = condition, fill = condition), alpha = 0.1) +
    scale_color_manual(values = c(ga_midgreen, ga_navy, ga_darkgreen)) +
    scale_fill_manual(values = c(ga_midgreen, ga_navy, ga_darkgreen)) +
    theme_bw() +
    ggtitle(paste0("Peak activity for ", peak, "\nNormalisation: ", normalisation))
}

#FILTER

#calculate mean values to use for filtering
means <- dm2 %>%
  group_by(condition) %>%
  summarise(mean = mean(value))

mean_16h <- means %>% filter(condition == col_16h) %>% pull(mean)
mean_5d <- means %>% filter(condition == col_5d) %>% pull(mean)

#require that Th0 status is higher than unstimulated
#AND require that at least one of the Th0 peaks has a normalised value of at least the mean
dm_filt <- dm %>%
  filter(!!sym(col_16h) > mean_16h | !!sym(col_5d) > mean_5d) %>%
  filter(!!sym(col_uns) < !!sym(col_16h) | !!sym(col_uns) < !!sym(col_5d)) %>%
  select(chr, start, end, !!sym(col_16h), !!sym(col_5d), !!sym(col_uns))

#SAVE
file_filt <- paste0(dir, peak, "_counts_", normalisation, "_filtered.txt") 
fwrite(dm_filt, file_filt, col.names = T, row.names = F, sep = "\t")

#EXTEND
#Extend the ATAC windows by 1, 3 or 5kb to make the overlap with the SNPs bigger

if(peak == "ATAC"){
  windows <- c(1000, 3000, 5000)
  
  for (window in windows){
    dm_ext <- dm_filt %>%
      mutate(start = start - window,
             end = end + window)
    
    print(summary(dm_ext$end-dm_ext$start))
    
    #SAVE
    file_filt_ext <- paste0(dir, peak, "_counts_", normalisation, "_extended_", window, "_filtered.txt") 
    fwrite(dm_ext, file_filt_ext, col.names = T, row.names = F, sep = "\t")
  }
}





