#' Visualize a matched design as a
#'
#' Visualize which units are connected to which other units after some stratification.


#' @param strata_indicator records membership in a strata. For example, it
#' could be the factor that is created by fullmatch or pairmatch from the
#' optmatch package. It should have names. The nodes will be labeled with these
#' names

#' @param node_color is a vector that will be used to color the nodes in the
#' graph. For example, it could be the treatment variable in a bipartite
#' matched design. It should be in the same order as the strata_indicator variable.

#' @return a ggraph object

#' @example inst/examples/stratification_graph.R

#' @import ggraph tidygraph
#' @export

stratification_graph <- function(strata_indicator, node_color) {
  require(ggraph)
  require(tidygraph)

  ## first make an adjacency matrix using the strata indicator
  s <- droplevels(strata_indicator[!is.na(strata_indicator)])
  z <- node_color[!is.na(strata_indicator)]
  names(z) <- names(s)

  adj_mat <- outer(s, s, FUN = function(x, y) {
    as.numeric(x == y)
  })

  graph_obj <- as_tbl_graph(adj_mat, directed = FALSE)
  graph_obj <- graph_obj %>%
    activate(nodes) %>%
    mutate(trt = z)

  ## I like fr and graphopt, but using "nicely" for now
  # layout0 <- create_layout(graph_obj0, layout = 'igraph', algorithm = 'fr')
  thelayout <- create_layout(graph_obj, layout = "igraph", algorithm = "nicely")
  g <- ggraph(thelayout) +
    geom_node_point(aes(fill = as.factor(trt)), show.legend = FALSE) +
    geom_edge_link2() +
    # geom_node_label(aes(label=name,color=trt)) +
    # geom_edge_diagonal(colour = "black") +
    geom_node_label(aes(label = name, colour = trt),
      repel = FALSE, show.legend = FALSE, label.r = unit(0.5, "lines"),
      label.padding = unit(.01, "lines"), label.size = 0
    ) +
    theme(legend.position = "none") +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid,
      panel.border = element_rect(fill = "transparent", color = "black"),
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )

  return(g)
}
