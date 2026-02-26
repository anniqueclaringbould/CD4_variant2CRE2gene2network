## Count assay signal at elements from bigwig file

# save.image("RDA/count_chromatin_assay_bigwig.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(GenomicRanges)
  library(rtracklayer)
})

# load element from bed file
elements <- read_tsv(snakemake@input$elements, col_names = c("chr", "start", "end", "name"),
                     show_col_types = FALSE)

# load assay signal from bigwig file
assay <- import(snakemake@input$assay)

# make GRanges object from elements
elements_gr <- makeGRangesFromDataFrame(elements, starts.in.df.are.0based = TRUE)
names(elements_gr) <- elements$name

# compute total signal in bigwig file
total_signal <- sum(mcols(assay)$score * width(assay))

# overlap elements with assay signal
overlaps <- findOverlaps(elements_gr, assay)
overlaps <- split(overlaps, from(overlaps))
names(overlaps) <- names(elements_gr[as.numeric(names(overlaps))])

# function to compute signal for one enhancer
count_signal <- function(overlap, elements, assay) {
  
  # get signal only overlapping element
  element <- elements[unique(queryHits(overlap))]
  signal <- assay[unique(subjectHits(overlap))]
  signal_element <- pintersect(signal, element)
  
  # sum up signal across element
  signal_count <- sum(signal_element$score * width(signal_element))
  
  return(signal_count)
  
}

# calculate assay signal at each element
signal <- vapply(overlaps, FUN.VALUE = numeric(1), FUN = count_signal, elements = elements_gr,
                 assay = assay)

# convert to table, add to elements and add total signal
elements <- signal %>% 
  enframe(name = "name", value = "signal_count") %>% 
  left_join(elements, ., by = "name") %>% 
  mutate(signal_total = total_signal)

# save output to file
write_tsv(elements, snakemake@output[[1]])
