library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
library(scbp) # install_github("rnabioco/scbp")

theme_set(theme_cowplot())

# Calculate Force Directed Layout ------------------------------------------------

#' Calculate force directed graph coordinates
#'
#' @param seurat_obj seurat_obj
#' @param k_neighbors number of nearest neighbors to compute
#' @param seed integer seed for reproducible analysis
#'
#' @return data.frame with columns necessary for plotting with ggnetwork,
#' and all columns from meta.data slot
#'
#' @importFrom ggnetwork ggnetwork
#' @importFrom RANN nn2
#' @export
calc_graph <- function(sobj, k_neighbors = 15, seed = 42) {

  mat <- sobj@reductions$pca@cell.embeddings

  # calculate nearest neighbors in PCA space
  knn.info <- RANN::nn2(mat, k = k_neighbors)
  knn <- knn.info$nn.idx

  # convert to adjacency matrix
  adj <- matrix(0, nrow(mat), nrow(mat))
  rownames(adj) <- colnames(adj) <- rownames(mat)
  for(i in seq_len(nrow(mat))) {
    adj[i,rownames(mat)[knn[i,]]] <- 1
  }

  # make force directed graph directly from adjacency matrix
  set.seed(seed)
  gn <- ggnetwork::ggnetwork(adj,
                             layout = "fruchtermanreingold",
                             niter = 500)

  mdata <- tibble::as_tibble(sobj@meta.data, rownames = "cell")

  # add in metadata
  gn$vertex.names <- as.character(gn$vertex.names)
  gn <- left_join(gn, mdata, by = c("vertex.names" = "cell"))

  gn

}

#' Plot force directed graph with cell annotations
#'
#' @param graph_df data.frame produced by calc_graph
#' @param color_by column from graph_df to color cells by
#'
#' @return ggplot object
#'
#' @import ggnetwork
#' @export
plot_graph <- function(graph_df,
                       color_by = "orig.ident"){

  graph_df <- graph_df %>% arrange_at(.vars = color_by)

  p <- ggplot(graph_df,
              aes(x, y, xend = xend, yend = yend)) +
    geom_edges(
      aes_string(color = color_by),
      alpha = 0.1,
      arrow = arrow(length = unit(0.1, "pt"),
                    type = "closed"),
      curvature = 0.05) +
    geom_nodes(aes_string(color = color_by),
               size = 0.1)

  if(is_discrete(graph_df[[color_by]])){
    p <- p + scale_color_manual(color_by, values = discrete_palette_default)
  } else {
    max_y <- c(0, max(graph_df[[color_by]]))
    cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])

    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = cols,
                                   name = color_by)
  }

  p +
    theme_blank() +
    theme(legend.position = "bottom")
}



is_discrete <- function(x){
  if (typeof(x) %in% c(
    "character",
    "logical"
  ) | is.factor(x)) {
    discrete <- TRUE
  } else {
    discrete <- FALSE
  }
  discrete
}



