##############################################
## Script to combine significant DE results ##
######## from NT calibration results #########
##############################################

#### SETUP ####
#libraries
library(dplyr)
library(tidyverse)
library(limma)
library(edgeR)
library(ggpubr)
library(Seurat)
library(BPCells)
library(Matrix)
library(data.table)
library(ggrastr)

#directories
root <- "/g/steinmetz/project/otar_2063/promoter_screen/"
indir <- paste0(root, "interpretation/DE/output/per_gene/non_targeting/")
outdir <- paste0(root, "interpretation/DE/output/combined/non_targeting/")

#### READ DATA ####
file_name <- "DEG_single_cell_limma_fakePerturbation_[0-9]+_vs_all"
files <- list.files(path = indir, pattern = file_name)

res <- rbindlist(
  lapply(files, function(file){
    
    #extract fake perturbation number from filename
    pert <- sub(".*fakePerturbation_([0-9]+)_vs_all_.*", "\\1", file)
    
    #read file
    dm <- fread(paste0(indir, file)) %>%
      mutate(Perturbation = pert) %>%
      dplyr::select(Perturbation, Target = "V1", everything())
    
  }),
  use.names = TRUE, fill = TRUE
)

#### ADJUST ####

res <- res %>%
  mutate(significant = ifelse(adj.P.Val < 0.1, TRUE, FALSE)) %>%
  distinct()

filename <- gsub("_\\[0-9\\]\\+", "", file_name)
fwrite(res, paste0(outdir, filename, "_combined_p01_", Sys.Date(), ".tsv.gz"))

res_sig <- res %>%
  filter(significant == TRUE) 

res_sig %>%
  group_by(Target) %>%
  tally() %>%
  arrange(-n)

res_sig %>%
  group_by(Perturbation) %>%
  tally() %>%
  arrange(-n)

#### PLOT ####

#Volcano plots
p1 <- res %>%
  filter(P.Value < 0.001) %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), col = significant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

p2 <- res_sig %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), col = AveExpr)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

pdf(paste0(outdir, filename, "_Volcano_plots.pdf"))
rasterize(p1, dpi = 300)
rasterize(p2, dpi = 300)
dev.off()

#QQ-plots
#using code from https://slowkow.com/notes/ggplot2-qqplot/

#calculate inflation
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

inflation(res$P.Value)

#percentage of positive results
nrow(res_sig) / nrow(res) * 100

#number of positive results per 100.000 tests
nrow(res_sig) / nrow(res) * 100000

#QQ-plot
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

qq <- gg_qqplot(res$P.Value) +
  theme_classic()
qq.r <- rasterise(qq, dpi = 700)
ggsave(plot = qq.r, paste0(outdir, filename, "_QQ-plot.pdf"))

