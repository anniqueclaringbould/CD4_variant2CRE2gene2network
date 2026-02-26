## Make sceptre QQ plots for all CRE screen panels from sceptre objects with stored output

# save.image("RDA/make_cre_screen_qq_plots.rda")
# stop()

suppressPackageStartupMessages({
  library(sceptre)
  library(tidyverse)
  library(cowplot)
})

# Define functions ---------------------------------------------------------------------------------

# function to load one sceptre object and make discovery analysis plots
plot_discovery_analysis <- function(sceptre_object, name) {
  
  message("Making plots for: ", sceptre_object)
  
  # load sceptre object
  sceptre_object <- readRDS(sceptre_object)
  
  # make individual discovery analysis plot as ggplot objects that can be further altered
  disc_plots <- plot(sceptre_object, return_indiv_plots = TRUE)
  
  # manually set new titles for each plot to include the specified name associated with the data
  disc_plots <- lapply(disc_plots, FUN = function(x) x + labs(title = name))
  
  # free unused memory (defensive, probably not needed)
  rm(sceptre_object)
  invisible(gc())
  
  return(disc_plots)
  
}

# arrange plots into one figure
arrange_figure <- function(plots, nrow = NULL, ncol = NULL, title, title_size = 10, scale = 1, 
                           legend = NULL) {
  
  # create title object
  title_plot <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0, 
      size = title_size
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  # arrange plots and add title
  plots <- plot_grid(plotlist = plots, nrow = nrow, ncol = ncol, scale = scale)
  plots <- plot_grid(title_plot, plots, ncol = 1, rel_heights = c(0.1, 1))
  
  # add legend if provided
  if (!is.null(legend)) {
    plots <- plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.1))
  }
  
  return(plots)
  
}


# Make plots ---------------------------------------------------------------------------------------

# sceptre objects for all panels
sceptre_objects <- unlist(snakemake@config$finished_sceptre_objects)
names(sceptre_objects) <- sub("Panel", "Panel ", names(sceptre_objects))

# order sceptre objects numerically
sceptre_objects <- sceptre_objects[order(as.numeric(sub("Panel ", "", names(sceptre_objects))))]

# make discovery analysis plots for all sceptre objects
disc_plots <- lapply(names(sceptre_objects), FUN = function(x) {
  plot_discovery_analysis(sceptre_objects[[x]], name = x)
})

# get plots across panels for each plot type
qq_plots_bulk <- lapply(disc_plots, FUN = "[[", 1)
qq_plots_tail <- lapply(disc_plots, FUN = "[[", 2)

# extract legend from bulk plot
bulk_legend <- get_legend(qq_plots_bulk[[1]])

## Make figure with 2x8 layout ---------------------------------------------------------------------

# manually set axis titles for qq_plots_tail plots
qq_plots_tail_28 <- qq_plots_tail
qq_plots_tail_28[c(1, 9)] <- lapply(qq_plots_tail_28[c(1, 9)], FUN = function(p) {
  p + theme(axis.title.y = element_text(angle = 90))
})
qq_plots_tail_28[1:8] <- lapply(qq_plots_tail_28[1:8], FUN = function(p) {
  p + theme(axis.title.x = element_blank())
})

# combine tail QQ plots into one multi-panel figure (8 x 2)
qq_plots_tail_28 <- arrange_figure(qq_plots_tail_28, ncol = 8, title = "QQ plots (tail)",
                                   legend = bulk_legend)

# save plots to files
ggsave(qq_plots_tail_28, file = snakemake@output$qq_plots_tail_2by8, height = 6, width = 20)

## Make figure with 4x4 layout ---------------------------------------------------------------------

# manually set axis titles for qq_plots_tail plots
qq_plots_tail_44 <- qq_plots_tail
qq_plots_tail_44[c(1, 5, 9, 13)] <- lapply(qq_plots_tail_44[c(1, 5, 9, 13)], FUN = function(p) {
  p + theme(axis.title.y = element_text(angle = 90))
})
qq_plots_tail_44[1:12] <- lapply(qq_plots_tail_44[1:12], FUN = function(p) {
  p + theme(axis.title.x = element_blank())
})

# combine tail QQ plots into one multi-panel figure (8 x 2)
qq_plots_tail_44 <- arrange_figure(qq_plots_tail_44, ncol = 4, title = "QQ plots (tail)",
                                   legend = bulk_legend)

# save plots to files
ggsave(qq_plots_tail_44, file = snakemake@output$qq_plots_tail_4by4, height = 10, width = 10)
