#' Visualize a matched design as graph
#'
#' Visualize which units are connected to which other units after some
#' stratification using nodes and edges.

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



#' Visualize the distances that remain after creating a stratified design
#
#' Visualize distances between units within sets after matching using a violin
#' plot, a boxplot, and dots. This function works best with covariates that
#' have more than two values. The sets with the largest and smallest absolute
#' differences are labeled by default.

#' @param strata_indicator records membership in a strata. For example, it
#' could be the factor that is created by fullmatch or pairmatch from the
#' optmatch package. It should have names.

#' @param trt is a vector that indicates he treatment variable in a bipartite
#' matched design. It should be in the same order as the strata_indicator variable.

#' @param covariate is a vector with the same names as the strata_indicator and
#' node_color. The function will display distances within sets on this covariate.

#' @return a ggplot object

#' @example inst/examples/covariate_matched_dist_plot.R

#' @import ggplot2 dplyr
#' @export

covariate_matched_dist_plot <- function(strata_indicator, trt, covariate) {
  stopifnot(all.equal(length(strata_indicator), length(covariate)))
  stopifnot(all.equal(length(strata_indicator), length(trt)))

  dat <- na.omit(data.frame(s = strata_indicator, x = covariate, z = trt))

  cov_diffs_dat <- dat %>%
    group_by(s) %>%
    reframe(abs_cov_diffs = abs(x[z == 1] - x[z == 0])) %>%
    ungroup()

  cov_diffs_to_label <- cov_diffs_dat %>% filter(abs_cov_diffs %in% range(abs_cov_diffs))

  data_summary <- function(x) {
    ptiles <- quantile(x, c(.1, .5, .9))
    return(c(y = ptiles[["50%"]], ymin = ptiles[["10%"]], ymax = ptiles[["90%"]]))
  }

  ## For now, display for only one design so x is just a constant and we don't
  ## display the x-axis

  g <- ggplot(data = cov_diffs_dat, aes(
    y = abs_cov_diffs,
    x = rep(1, nrow(cov_diffs_dat))
  )) +
    geom_violin(aes(fill = NULL)) +
    geom_boxplot(width = .1, aes(fill = NULL)) +
    stat_summary(fun.data = data_summary, aes(fill = NULL)) +
    ylab("Within Set Abs Covariate Differences") +
    geom_dotplot(
      # fill = as.numeric(cov_diffs_dat$s),
      binaxis = "y", stackdir = "center", stackratio = .7,
      dotsize = .8, method = "histodot", binwidth = 1
    ) +
    geom_text(
      data = cov_diffs_to_label,
      aes(y = abs_cov_diffs, x = rep(1, nrow(cov_diffs_to_label)), label = s),
      check_overlap = TRUE,
      fontface = 2,
      nudge_x = .25
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(g)
}
