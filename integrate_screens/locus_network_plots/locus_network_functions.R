## Functions to create locus network plots
## Includes optional edge bundling between perturb_seq layers
## Bundles edges between all consecutive perturb_seq layers to avoid overplotting
## Intra-layer edges curve outward for better visibility

suppressPackageStartupMessages({
  library(igraph)
})

#' Create a locus network plot
#'
#' Generates a hierarchical network visualization showing CRE-gene interactions from TAP-seq data
#' and downstream gene-gene interactions from Perturb-seq data. The network is displayed as a tree
#' with CREs at the top, their target genes in the second layer, and downstream regulatory effects
#' in subsequent layers.
#'
#' @param enh_hits A data frame containing CRE TAP-seq screen results with columns including
#'   'regulator' (CRE identifier), 'target' (gene name), 'pert_chr', 'pert_start', 'pert_end'
#'   (CRE coordinates), and 'effect_size'.
#' @param prom_hits A data frame containing promoter Perturb-seq screen results with columns
#'   including 'regulator', 'target', and 'effect_size'.
#' @param start_gene Character vector of gene name(s) to use as the starting point for building
#'   the network. These should be target genes of CREs in enh_hits.
#' @param layers Integer specifying the number of downstream Perturb-seq layers to include in the
#'   network. Default is 2.
#' @param directed Logical indicating whether to create a directed graph. Default is TRUE.
#' @param jitter_amount Numeric value controlling the amount of vertical jittering applied to nodes
#'   in Perturb-seq layers to reduce overlap. Default is 0.3.
#' @param highlight_genes Character vector of gene names to highlight in the plot. Highlighted genes
#'   are displayed with larger nodes and darker color. Default is NULL.
#' @param return_graph Logical indicating whether to return the igraph object instead of just
#'   plotting. When TRUE, returns the graph with layout coordinates stored as vertex attributes.
#'   Default is FALSE.
#' @param node_size Numeric value specifying the base size of nodes. Default is 5.
#' @param cluster_gap Numeric value controlling the horizontal gap between clusters of nodes that
#'   share different parent regulators in the most downstream layer. Default is 0.5.
#' @param min_node_spacing Numeric value specifying the minimum horizontal spacing between nodes
#'   within a layer. Default is 0.25.
#' @param legend_cex Numeric value controlling the font size of the gene legend below the plot.
#'   Default is 0.5.
#' @param bundle_edges Logical indicating whether to bundle edges between consecutive Perturb-seq
#'   layers. Bundling groups edges from the same source node together to reduce visual clutter.
#'   Default is TRUE.
#' @param bundle_tension Numeric value between 0 and 1 controlling how tightly bundled edges
#'   converge. Higher values create tighter bundles. Default is 0.5.
#' @param intra_layer_curve Numeric value controlling the curvature of edges between nodes within
#'   the same layer. Default is 0.3.
#'
#' @return If return_graph is FALSE (default), the function plots the network and returns NULL
#'   invisibly. If return_graph is TRUE, returns an igraph object with vertex attributes
#'   'X Location', 'Y Location', 'Fill Color', and 'layer_name' for use in downstream analysis
#'   or export.
#'
#' @details
#' The network layout is organized hierarchically:
#' \itemize{
#'   \item Top layer: CREs (cis-regulatory elements), colored red
#'   \item Second layer: CRE target genes, colored light blue
#'   \item Lower layers: Downstream genes from Perturb-seq, numbered and shown in a legend
#' }
#'
#' Nodes are ordered within layers based on genomic coordinates (CREs and genes) or parent node
#' positions (downstream layers). Edge colors indicate effect direction: teal for negative effects
#' and orange for positive effects.
#'
#' @examples
#' \dontrun{
#' # Basic network plot
#' make_locus_network(enh_hits, prom_hits, start_gene = "IL2RA")
#'
#' # Network with 2 downstream layers and highlighted genes
#' make_locus_network(enh_hits, prom_hits, start_gene = "IL2RA", layers = 2,
#'                    highlight_genes = c("FOXP3", "CTLA4"))
#'
#' # Return graph object for export
#' g <- make_locus_network(enh_hits, prom_hits, start_gene = "IL2RA", return_graph = TRUE)
#' }
#'
#' @export
make_locus_network <- function(enh_hits, prom_hits, start_gene, layers = 2, directed = TRUE,
                               jitter_amount = 0.3, highlight_genes = NULL, return_graph = FALSE,
                               node_size = 5, cluster_gap = 0.5, min_node_spacing = 0.25,
                               legend_cex = 0.5, bundle_edges = TRUE, bundle_tension = 0.5,
                               intra_layer_curve = 0.3) {

  # get significant enhancers for given genes
  enh_hits_gene <- filter(enh_hits, target %in% start_gene)

  # get perturb-seq hits across the specified number of layers for the given gene
  prom_hit_layers <- get_perturb_seq_hits(prom_hits, start_gene = enh_hits_gene$target,
                                          layers = layers)

  # combine promoter and enhancer hits
  hits <- bind_rows(enh_hits_gene, prom_hit_layers)

  # create graph
  graph <- graph_from_data_frame(hits, directed = directed)

  # create colors for nodes
  node_colors <- rep("#bfc8e3", vcount(graph))
  names(node_colors) <- V(graph)$name
  node_colors[start_gene] <- "#bfc8e3"
  node_colors[enh_hits_gene$regulator] <- "#e57473"

  # create sizes for nodes
  node_sizes <- rep(node_size, vcount(graph))
  names(node_sizes) <- V(graph)$name
  node_label_cex <- unique(node_sizes) * 0.15

  # create distances for labels from nodes
  node_dists <- rep(node_size * 0.1, vcount(graph))
  names(node_dists) <- V(graph)$name

  # highlight specified genes in red with larger size
  if (!is.null(highlight_genes)) {
    highlight_genes <- intersect(highlight_genes, names(node_colors))
    node_colors[highlight_genes] <- "#334068"
    node_sizes[highlight_genes] <- unique(node_sizes) * 1.6
    node_dists[highlight_genes] <- unique(node_dists) * 1.6
  }

  # create layout
  layout <- layout_as_tree(graph, root = enh_hits_gene$regulator)
  node_names <- V(graph)$name
  original_y <- layout[, 2]  # store original layer assignments used later

  # get TAP-seq layer nodes
  tap_seq_layers <- sort(unique(layout[, 2]), decreasing = TRUE)[c(1,2)]

  # order CRE nodes (top layer) by perturbation coordinates
  cre_layer <- tap_seq_layers[1]
  cre_indices <- which(layout[, 2] == cre_layer)

  if (length(cre_indices) > 1) {
    cre_node_names <- node_names[cre_indices]

    # get CRE coordinates from hits table
    cre_coords <- hits %>%
      filter(regulator %in% cre_node_names) %>%
      select(regulator, pert_chr, pert_start, pert_end) %>%
      distinct()

    # only reorder if we have valid coordinates
    if (nrow(cre_coords) > 0 && any(!is.na(cre_coords$pert_chr))) {
      cre_coords <- cre_coords %>%
        filter(!is.na(pert_chr)) %>%
        arrange(pert_chr, pert_start, pert_end)

      sorted_x <- sort(layout[cre_indices, 1])

      for (i in seq_len(nrow(cre_coords))) {
        node_idx <- which(node_names == cre_coords$regulator[i])
        layout[node_idx, 1] <- sorted_x[i]
      }
    }
  }

  # order gene nodes (second layer) by gene TSS coordinates
  gene_layer <- tap_seq_layers[2]
  gene_indices <- which(layout[, 2] == gene_layer)

  if (length(gene_indices) > 1) {
    gene_node_names <- node_names[gene_indices]

    # get gene coordinates from hits table
    gene_coords <- hits %>%
      filter(target %in% gene_node_names) %>%
      select(target, gene_chr, gene_tss) %>%
      distinct()

    # only reorder if we have valid coordinates
    if (nrow(gene_coords) > 0 && any(!is.na(gene_coords$gene_chr))) {
      gene_coords <- gene_coords %>%
        filter(!is.na(gene_chr)) %>%
        arrange(gene_chr, gene_tss)

      sorted_x <- sort(layout[gene_indices, 1])

      for (i in seq_len(nrow(gene_coords))) {
        node_idx <- which(node_names == gene_coords$target[i])
        layout[node_idx, 1] <- sorted_x[i]
      }
    }
  }

  # order nodes in perturb-seq layers based on parent node x-coordinates
  perturb_seq_layers <- sort(unique(layout[, 2]), decreasing = TRUE)[-c(1,2)]

  # identify the most downstream layer (smallest y value)
  if (length(perturb_seq_layers) > 0) {
    downstream_layer <- min(perturb_seq_layers)
  } else {
    downstream_layer <- NULL
  }

  for (layer_y in perturb_seq_layers) {
    layer_indices <- which(layout[, 2] == layer_y)

    if (length(layer_indices) > 1) {
      layer_node_names <- node_names[layer_indices]

      # calculate mean x-coordinate of parent nodes for each node in this layer
      parent_x_means <- sapply(layer_node_names, function(node) {
        # get parent nodes (regulators) from the graph
        parents <- neighbors(graph, node, mode = "in")
        if (length(parents) > 0) {
          parent_indices <- which(node_names %in% names(parents))
          mean(layout[parent_indices, 1])
        } else {
          layout[which(node_names == node), 1]
        }
      })

      # order nodes by parent x-coordinate means
      node_order <- order(parent_x_means)
      sorted_x <- sort(layout[layer_indices, 1])

      # assign new x positions based on parent ordering
      for (i in seq_along(node_order)) {
        node_idx <- which(node_names == layer_node_names[node_order[i]])
        layout[node_idx, 1] <- sorted_x[i]
      }

      # for the most downstream layer, add gaps between clusters from different parents
      if (!is.null(downstream_layer) && layer_y == downstream_layer && cluster_gap > 0) {
        layout <- add_cluster_gaps(graph, layout, node_names, layer_indices,
                                   layer_node_names, cluster_gap, min_node_spacing)
      }
    }
  }

  # determine jittering strategy based on bundling
  if (!identical(unique(layout[, 2]), layout[, 2])) {

    if (bundle_edges && length(perturb_seq_layers) >= 2 && jitter_amount > 0) {
      # when bundling is enabled, jitter clusters together (not individual nodes)
      # apply cluster jittering to all perturb_seq layers
      for (layer_y in perturb_seq_layers) {
        layout <- jitter_clusters(graph, layout, node_names, layer_y, jitter_amount)
      }

    } else if (length(perturb_seq_layers) > 1) {
      # no bundling: jitter all layers after layer 1
      layers_to_jitter <- perturb_seq_layers[-1]
      layout <- jitter_layers(layout, layers = layers_to_jitter, jitter_amount = jitter_amount)
    }

    # normalize layout to [-1, 1] (same as igraph default)
    layout_norm <- layout
    layout_norm[, 1] <- 2 * (layout[, 1] - min(layout[, 1])) /
      (max(layout[, 1]) - min(layout[, 1])) - 1
    layout_norm[, 2] <- 2 * (layout[, 2] - min(layout[, 2])) /
      (max(layout[, 2]) - min(layout[, 2])) - 1

  } else {
    layout_norm <- layout
  }

  # add layer names as node attributes to graph
  V(graph)$layer_name <- case_when(
    original_y == tap_seq_layers[[1]] ~ "CRE",
    original_y == tap_seq_layers[[2]] ~ "CRE-gene",
    original_y %in% perturb_seq_layers ~ paste0("Hop", abs(original_y - max(original_y)) - 1)
  )

  # create node labels - use numbers inside nodes for perturb_seq_layer nodes
  node_labels <- V(graph)$name
  perturb_seq_indices <- which(original_y %in% perturb_seq_layers)
  perturb_seq_node_names <- node_names[perturb_seq_indices]

  # create number labels for perturb_seq nodes and store mapping for legend
  if (length(perturb_seq_indices) > 0) {
    # order by x-coordinate for consistent numbering
    x_order <- order(layout_norm[perturb_seq_indices, 1])
    ordered_names <- perturb_seq_node_names[x_order]

    # create mapping: number -> gene name (exclude highlight genes)
    non_highlight_names <- setdiff(ordered_names, highlight_genes)
    non_highlight_numbers <- which(ordered_names %in% non_highlight_names)
    legend_mapping <- data.frame(
      number = seq_along(non_highlight_numbers),
      gene = ordered_names[non_highlight_numbers],
      stringsAsFactors = FALSE
    )

    # assign number labels to perturb_seq nodes (except highlight genes)
    number_counter <- 1
    for (i in seq_along(ordered_names)) {
      node_idx <- which(node_names == ordered_names[i])

      if (!is.null(highlight_genes) && ordered_names[i] %in% highlight_genes) {
        # highlight genes: keep gene name label above node
        node_labels[node_idx] <- ordered_names[i]
        node_dists[node_idx] <- 0.8  # label outside node
      } else {
        # regular nodes: number inside node
        node_labels[node_idx] <- as.character(number_counter)
        node_dists[node_idx] <- 0  # center label inside node
        number_counter <- number_counter + 1
      }
    }
  } else {
    legend_mapping <- NULL
  }

  # create edge colors based on effect direction
  edge_colors <- ifelse(E(graph)$effect_size < 0, "#07bab9", "#e78b24")

  # get edge list for identifying special edges

  edge_list <- as_edgelist(graph)

  # identify edges between consecutive perturb_seq layers for bundling
  # (layer 1 to 2, layer 2 to 3, etc.)
  bundled_edge_indices <- c()
  if (bundle_edges && length(perturb_seq_layers) >= 2) {
    for (layer_idx in seq_len(length(perturb_seq_layers) - 1)) {
      source_layer <- perturb_seq_layers[layer_idx]
      target_layer <- perturb_seq_layers[layer_idx + 1]
      source_layer_nodes <- node_names[which(original_y == source_layer)]
      target_layer_nodes <- node_names[which(original_y == target_layer)]

      for (i in seq_len(nrow(edge_list))) {
        from_node <- edge_list[i, 1]
        to_node <- edge_list[i, 2]
        if (from_node %in% source_layer_nodes && to_node %in% target_layer_nodes) {
          bundled_edge_indices <- c(bundled_edge_indices, i)
        }
      }
    }
    bundled_edge_indices <- unique(bundled_edge_indices)
  }

  # identify intra-layer edges (edges within the same layer) - only when bundling
  intra_layer_edge_indices <- c()
  if (bundle_edges) {
    for (i in seq_len(nrow(edge_list))) {
      from_node <- edge_list[i, 1]
      to_node <- edge_list[i, 2]
      from_idx <- which(node_names == from_node)
      to_idx <- which(node_names == to_node)
      # check if both nodes are in the same original layer
      if (length(from_idx) > 0 && length(to_idx) > 0 &&
          original_y[from_idx] == original_y[to_idx]) {
        intra_layer_edge_indices <- c(intra_layer_edge_indices, i)
      }
    }
  }

  # combine all edges to draw manually (bundled + intra-layer)
  manual_edge_indices <- unique(c(bundled_edge_indices, intra_layer_edge_indices))

  # plot graph - suppress manually drawn edges
  if (length(manual_edge_indices) > 0) {
    # create a modified graph without the manual edges for initial plot
    graph_no_manual <- delete_edges(graph, manual_edge_indices)
    edge_colors_no_manual <- edge_colors[-manual_edge_indices]

    plot(graph_no_manual,
         layout = layout_norm[match(V(graph_no_manual)$name, node_names), ],
         rescale = FALSE,
         vertex.size = node_sizes[V(graph_no_manual)$name],
         vertex.label = node_labels[match(V(graph_no_manual)$name, node_names)],
         vertex.label.color = "black",
         vertex.label.cex = node_label_cex,
         vertex.label.dist = node_dists[V(graph_no_manual)$name],
         vertex.label.degree = -pi/2,
         vertex.color = node_colors[V(graph_no_manual)$name],
         edge.color = edge_colors_no_manual,
         edge.arrow.size = 0.5)

    # draw bundled edges manually
    if (bundle_edges && length(bundled_edge_indices) > 0) {
      draw_bundled_edges(graph, layout_norm, node_names, bundled_edge_indices,
                         edge_colors, bundle_tension, arrow_size = 0.5)
    }

    # draw intra-layer edges with outward curves
    if (length(intra_layer_edge_indices) > 0) {
      draw_intra_layer_edges(graph, layout_norm, node_names, original_y,
                             intra_layer_edge_indices, edge_colors,
                             curve_amount = intra_layer_curve, arrow_size = 0.5)
    }
  } else {
    plot(graph,
         layout = layout_norm,
         rescale = FALSE,
         vertex.size = node_sizes,
         vertex.label = node_labels,
         vertex.label.color = "black",
         vertex.label.cex = node_label_cex,
         vertex.label.dist = node_dists,
         vertex.label.degree = -pi/2,
         vertex.color = node_colors,
         edge.color = edge_colors,
         edge.arrow.size = 0.5)
  }

  # add box around nodes of perturb-seq layers
  for (i in seq_along(perturb_seq_layers)) {
    draw_layer_box(layout_norm, layer = perturb_seq_layers[[i]], original_y = original_y,
                   layer_title = i, x_min = -1, x_max = 1, padding = 0.1)
  }

  # add legend for numbered nodes below the plot
  if (!is.null(legend_mapping) && nrow(legend_mapping) > 0) {
    draw_node_legend(legend_mapping, y_start = -1.22, x_start = -1.15, x_end = 1.15, n_cols = 8,
                     text_cex = legend_cex, row_height = 0.05)
  }

  # save network graph if specified
  if (return_graph == TRUE) {
    V(graph)$`X Location` <- layout_norm[, 1] * 500
    V(graph)$`Y Location` <- layout_norm[, 2] * -500
    V(graph)$`Fill Color` <- node_colors[V(graph)$name]
    return(graph)
  }

}


# function to draw bundled edges between layers
draw_bundled_edges <- function(graph, layout_norm, node_names, edge_indices, edge_colors,
                               tension = 0.5, arrow_size = 0.5) {

  edge_list <- as_edgelist(graph)

  # group edges by source node for bundling
  edge_groups <- list()
  for (i in edge_indices) {
    from_node <- edge_list[i, 1]
    if (is.null(edge_groups[[from_node]])) {
      edge_groups[[from_node]] <- list(indices = c(), to_nodes = c())
    }
    edge_groups[[from_node]]$indices <- c(edge_groups[[from_node]]$indices, i)
    edge_groups[[from_node]]$to_nodes <- c(edge_groups[[from_node]]$to_nodes, edge_list[i, 2])
  }

  # draw bundled edges for each source node
  for (from_node in names(edge_groups)) {
    group <- edge_groups[[from_node]]
    from_idx <- which(node_names == from_node)
    from_x <- layout_norm[from_idx, 1]
    from_y <- layout_norm[from_idx, 2]

    # calculate bundle point (where edges converge below source node)
    to_indices <- sapply(group$to_nodes, function(n) which(node_names == n))
    to_x_mean <- mean(layout_norm[to_indices, 1])
    to_y <- layout_norm[to_indices[1], 2]

    # bundle point is between source and targets
    bundle_y <- from_y + tension * (to_y - from_y)
    bundle_x <- from_x + tension * (to_x_mean - from_x)

    # draw each edge as a curve through the bundle point
    for (j in seq_along(group$indices)) {
      edge_idx <- group$indices[j]
      to_node <- group$to_nodes[j]
      to_idx <- which(node_names == to_node)
      to_x <- layout_norm[to_idx, 1]
      to_y_actual <- layout_norm[to_idx, 2]

      color <- edge_colors[edge_idx]

      # draw curved edge using bezier-like segments
      draw_bundled_edge_curve(from_x, from_y, bundle_x, bundle_y,
                              to_x, to_y_actual, color, arrow_size)
    }
  }
}


# function to draw a single bundled edge as a smooth curve
draw_bundled_edge_curve <- function(x0, y0, cx, cy, x1, y1, color, arrow_size = 0.5) {

  # create bezier-like curve using quadratic interpolation
  n_points <- 30
  t_vals <- seq(0, 1, length.out = n_points)

  # quadratic bezier curve: P(t) = (1-t)^2 * P0 + 2*(1-t)*t * C + t^2 * P1
  x_curve <- (1 - t_vals)^2 * x0 + 2 * (1 - t_vals) * t_vals * cx + t_vals^2 * x1
  y_curve <- (1 - t_vals)^2 * y0 + 2 * (1 - t_vals) * t_vals * cy + t_vals^2 * y1

  # draw the curve
  lines(x_curve, y_curve, col = color, lwd = 1)

  # draw arrowhead at the end
  arrow_length <- 0.05 * arrow_size
  # calculate angle at the end of the curve
  dx <- x_curve[n_points] - x_curve[n_points - 1]
  dy <- y_curve[n_points] - y_curve[n_points - 1]
  angle <- atan2(dy, dx)

  # arrowhead points
  arrow_angle <- pi / 6  # 30 degrees
  x_arrow1 <- x1 - arrow_length * cos(angle - arrow_angle)
  y_arrow1 <- y1 - arrow_length * sin(angle - arrow_angle)
  x_arrow2 <- x1 - arrow_length * cos(angle + arrow_angle)
  y_arrow2 <- y1 - arrow_length * sin(angle + arrow_angle)

  polygon(c(x1, x_arrow1, x_arrow2), c(y1, y_arrow1, y_arrow2), col = color, border = color)
}


# function to draw intra-layer edges with outward curves
draw_intra_layer_edges <- function(graph, layout_norm, node_names, original_y,
                                   edge_indices, edge_colors, curve_amount = 0.3,
                                   arrow_size = 0.5) {

  edge_list <- as_edgelist(graph)

  for (i in edge_indices) {
    from_node <- edge_list[i, 1]
    to_node <- edge_list[i, 2]

    from_idx <- which(node_names == from_node)
    to_idx <- which(node_names == to_node)

    from_x <- layout_norm[from_idx, 1]
    from_y <- layout_norm[from_idx, 2]
    to_x <- layout_norm[to_idx, 1]
    to_y <- layout_norm[to_idx, 2]

    color <- edge_colors[i]

    # determine curve direction based on layer position
    # curve downward (negative y direction) for most layers
    layer_y <- original_y[from_idx]
    all_layers <- sort(unique(original_y), decreasing = TRUE)

    # curve outward from the center of the plot (downward for upper layers)
    curve_direction <- -1  # default: curve downward

    # calculate control point for the curve (midpoint with vertical offset)
    mid_x <- (from_x + to_x) / 2
    mid_y <- (from_y + to_y) / 2

    # scale curve amount by horizontal distance for consistent appearance
    h_dist <- abs(to_x - from_x)
    scaled_curve <- curve_amount * max(h_dist, 0.2)

    # control point is offset perpendicular to the edge
    ctrl_y <- mid_y + curve_direction * scaled_curve

    # draw the curved edge
    draw_intra_layer_curve(from_x, from_y, mid_x, ctrl_y, to_x, to_y,
                           color, arrow_size)
  }
}


# function to draw a single intra-layer edge as a smooth curve
draw_intra_layer_curve <- function(x0, y0, cx, cy, x1, y1, color, arrow_size = 0.5) {

  # create bezier-like curve using quadratic interpolation
  n_points <- 30
  t_vals <- seq(0, 1, length.out = n_points)

  # quadratic bezier curve: P(t) = (1-t)^2 * P0 + 2*(1-t)*t * C + t^2 * P1
  x_curve <- (1 - t_vals)^2 * x0 + 2 * (1 - t_vals) * t_vals * cx + t_vals^2 * x1
  y_curve <- (1 - t_vals)^2 * y0 + 2 * (1 - t_vals) * t_vals * cy + t_vals^2 * y1

  # draw the curve
  lines(x_curve, y_curve, col = color, lwd = 1)

  # draw arrowhead at the end
  arrow_length <- 0.05 * arrow_size
  # calculate angle at the end of the curve
  dx <- x_curve[n_points] - x_curve[n_points - 1]
  dy <- y_curve[n_points] - y_curve[n_points - 1]
  angle <- atan2(dy, dx)

  # arrowhead points
  arrow_angle <- pi / 6  # 30 degrees
  x_arrow1 <- x1 - arrow_length * cos(angle - arrow_angle)
  y_arrow1 <- y1 - arrow_length * sin(angle - arrow_angle)
  x_arrow2 <- x1 - arrow_length * cos(angle + arrow_angle)
  y_arrow2 <- y1 - arrow_length * sin(angle + arrow_angle)

  polygon(c(x1, x_arrow1, x_arrow2), c(y1, y_arrow1, y_arrow2), col = color, border = color)
}


# function to add horizontal gaps between clusters of nodes from different parents
add_cluster_gaps <- function(graph, layout, node_names, layer_indices, layer_node_names,
                             cluster_gap, min_node_spacing) {

  # get parent info for each node (already ordered by x position)
  node_x_order <- order(layout[layer_indices, 1])
  ordered_nodes <- layer_node_names[node_x_order]

  parent_keys <- sapply(ordered_nodes, function(node) {
    parents <- neighbors(graph, node, mode = "in")
    paste(sort(names(parents)), collapse = "|")
  })

  # identify cluster boundaries using run-length encoding
  clusters <- rle(parent_keys)
  n_clusters <- length(clusters$lengths)

  # if only one cluster, no gaps needed
  if (n_clusters <= 1) return(layout)

  # get current x range
  current_x <- layout[layer_indices[node_x_order], 1]
  x_min <- min(current_x)
  x_max <- max(current_x)
  x_range <- x_max - x_min
  if (x_range == 0) x_range <- 2

  # distribute clusters evenly across the full width
  # each cluster gets space proportional to its size, plus gaps between clusters
  n_nodes <- length(ordered_nodes)
  n_gaps <- n_clusters - 1

  # calculate required width if using minimum spacing
  min_required_width <- (n_nodes - 1) * min_node_spacing

  # gap size between clusters
  gap_size <- cluster_gap * (x_range / n_clusters) * 1.5

  # if minimum required width plus gaps exceeds original range, expand the range
  total_gap_space <- n_gaps * gap_size
  if (min_required_width + total_gap_space > x_range) {
    x_range <- min_required_width + total_gap_space
  }

  node_space <- x_range - total_gap_space
  if (node_space < 0) node_space <- x_range * 0.5  # safety check

  # spacing between nodes within clusters (at least min_node_spacing)
  within_cluster_spacing <- max(node_space / max(n_nodes - 1, 1), min_node_spacing)

  # build new x positions: distribute clusters across full width with gaps
  new_x <- numeric(length(ordered_nodes))
  node_idx <- 1
  current_pos <- x_min

  for (cluster_idx in seq_len(n_clusters)) {
    cluster_size <- clusters$lengths[cluster_idx]
    cluster_end_idx <- node_idx + cluster_size - 1

    # place nodes within this cluster
    for (i in 0:(cluster_size - 1)) {
      new_x[node_idx + i] <- current_pos + i * within_cluster_spacing
    }

    # move position past this cluster
    current_pos <- current_pos + (cluster_size - 1) * within_cluster_spacing

    # add gap after this cluster (except for last cluster)
    if (cluster_idx < n_clusters) {
      current_pos <- current_pos + gap_size
    }

    node_idx <- cluster_end_idx + 1
  }

  # assign new positions back to layout
  for (i in seq_along(ordered_nodes)) {
    idx <- which(node_names == ordered_nodes[i])
    layout[idx, 1] <- new_x[i]
  }

  return(layout)
}


# function to get promoter hits for a given (set of) perturbations
get_promoter_hits <- function(prom_hits, perts, add_intra_set_hits = TRUE) {

  # get hits for all targets of perts
  hits_perts <- prom_hits %>%
    filter(regulator %in% unique(perts)) %>%
    filter(regulator != target)

  # add any hits between selected genes
  if (add_intra_set_hits == TRUE) {
    hits_perts <- prom_hits %>%
      filter(regulator %in% hits_perts$target) %>%
      filter(target %in% hits_perts$target) %>%
      filter(regulator != target) %>%
      bind_rows(hits_perts) %>%
      distinct()
  }

  return(hits_perts)

}

# function to get all perturbation targets perturb-seq layers
get_perturb_seq_hits <- function(prom_hits, start_gene, layers = 2) {

  # layers needs to be a positive integer and at least 1
  layers <- as.integer(layers)
  if (layers < 1) stop("At least one perturb-seq layer needs to be specified")

  # get promoter hits for initial gene
  prom_hits_gene <- list(get_promoter_hits(prom_hits, perts = start_gene))

  # get hits for subsequent layers
  for (i in seq_len(layers)[-1]) {
    prom_hits_layer <- get_promoter_hits(prom_hits, perts = prom_hits_gene[[i-1]]$target)
    if (nrow(prom_hits_layer) == 0) break
    prom_hits_gene[[i]] <- prom_hits_layer
  }

  # combine into one table and get unique interactions
  hits <- distinct(bind_rows(prom_hits_gene))

  return(hits)

}

# function to add vertical spacing within each layer by looping over layers
jitter_layers <- function(layout, layers, jitter_amount) {

  # add jittering to each layer specified by 'layers'
  for (layer in layers) {

    # get all nodes for given layer
    layer_nodes <- which(layout[, 2] == layer)
    n_nodes <- length(layer_nodes)

    # add y-coordinates jitter to nodes of that layer
    if (n_nodes > 1) {
      offsets <- seq(-jitter_amount, jitter_amount, length.out = n_nodes)
      layout[layer_nodes, 2] <- layout[layer_nodes, 2] + sample(offsets)
    }
  }

  return(layout)

}

# function to jitter clusters of nodes together (nodes sharing same parent get same offset)
jitter_clusters <- function(graph, layout, node_names, layer_y, jitter_amount) {

  # get nodes in the specified layer

  layer_indices <- which(layout[, 2] == layer_y)

  if (length(layer_indices) <= 1) return(layout)

  layer_node_names <- node_names[layer_indices]

  # group nodes by their parent(s)
  parent_keys <- sapply(layer_node_names, function(node) {
    parents <- neighbors(graph, node, mode = "in")
    if (length(parents) > 0) {
      paste(sort(names(parents)), collapse = "|")
    } else {
      node
    }
  })

  # get unique clusters
  unique_clusters <- unique(parent_keys)
  n_clusters <- length(unique_clusters)

  # if only one cluster, no jittering needed
  if (n_clusters <= 1) return(layout)

  # create offsets for each cluster (evenly spaced within jitter range)
  cluster_offsets <- seq(-jitter_amount, jitter_amount, length.out = n_clusters)

  # assign each cluster a jitter offset
  for (i in seq_along(unique_clusters)) {
    cluster_key <- unique_clusters[i]
    cluster_nodes <- layer_node_names[parent_keys == cluster_key]

    # apply the same offset to all nodes in this cluster
    for (node in cluster_nodes) {
      node_idx <- which(node_names == node)
      layout[node_idx, 2] <- layout[node_idx, 2] + cluster_offsets[i]
    }
  }

  return(layout)

}

# function to draw a box around a layer
draw_layer_box <- function(layout, layer, original_y, layer_title = NULL, x_min = -1, x_max = 1,
                           padding = 0.1) {

  # get all nodes for given layer
  layer_nodes <- which(original_y == layer)

  # title for box based on which layer it is
  box_title <- paste("Perturb-seq\nhop", layer_title)

  # draw box around nodes of given layer
  y_coords <- layout[layer_nodes, 2]
  y_mid <- (min(y_coords) + max(y_coords)) / 2
  rect(x_min - padding, min(y_coords) - padding,
       x_max + padding, max(y_coords) + padding,
       border = "gray40", lty = 2, lwd = 1.5)
  text(x_min - padding - 0.05, y_mid, box_title, adj = c(1, 0.5), cex = 0.7,
       col = "gray40")

}

# function to draw legend for numbered nodes
draw_node_legend <- function(legend_mapping, y_start = -1.22, x_start = -1.15, x_end = 1.15,
                             n_cols = 8, text_cex = 0.8, row_height = 0.05) {

  n_items <- nrow(legend_mapping)
  n_rows <- ceiling(n_items / n_cols)

  # calculate column width and row height
  col_width <- (x_end - x_start) / n_cols
  row_height <- row_height

  # draw legend title
  text(x_start, y_start + row_height, "Genes:", adj = c(0, 0.5), cex = text_cex * 1.2,
       col = "gray30", font = 2)

  # draw legend entries (fill by columns)
  for (i in seq_len(n_items)) {
    row <- ((i - 1) %% n_rows)
    col <- ((i - 1) %/% n_rows)

    x_pos <- x_start + col * col_width
    y_pos <- y_start - row * row_height

    label <- paste0(legend_mapping$number[i], ":", legend_mapping$gene[i])
    text(x_pos, y_pos, label, adj = c(0, 0.5), cex = text_cex, col = "gray30")
  }

}
