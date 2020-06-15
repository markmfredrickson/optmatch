##' @title Compute value of primal problem given flows and arc costs
##' @param distances An InfinitySparseMatrix giving distances
##' @param solution A MCFSolutions object
##' @return The value of the primal problem, i.e. sum of
##' products of \code{distances} with flow along arcs in \code{solution}
##' @author Hansen
##' @importFrom dplyr left_join
##' @importFrom dplyr filter 
evaluate_primal  <- function(distances, solution) {
    stopifnot(is(solution, "MCFSolutions"))
    anyflipped  <- any(solution@subproblems[["flipped"]])
    ## Following notation of Bertsekas *Network Optimization*, page 155,
    ## the primal problem value is
    ## \sum_{i,j} x_{ij} a_ij 
    ## where
    ##  - x_ij is the amount of flow along ij
    ##  - a_ij is the cost of the edge ij
    ## (Node prices don't have a role here.)

    nodes  <- as(nodeinfo(solution), "tbl_df")
    main_ij <- left_join(solution@arcs@matches,
                         dplyr::filter(nodes, upstream_not_down),
                         by = c("groups", "upstream" = "nodelabels")) %>%
               left_join(y = dplyr::filter(nodes, !upstream_not_down),
                         by = c("groups", "downstream" = "nodelabels"),
                         suffix = c(x = ".i", y = ".j")
                         )


    eld <- edgelist(distances)
    if (anyflipped)
        eld  <- rbind(eld,
                      edgelist(t(distances),
                               node.labels(solution)
                               )
                      )

    main_ij <- main_ij %>% left_join(eld,
                                     by = c("upstream" = "i",
                                            "downstream"= "j"),
                         suffix = c(x = "", y = ".dist")
                         )
    main_costs  <- sum(main_ij$dist)

    ## The below addresses a circumstance that doesn't
    ## currently arise, namely nonzero costs for arcs involving
    ## bookkeeping nodes.   It assumes that if/when that
    ## occurs:
    ## - bookkeeping node labels in the levels of the EdgeList's
    ##   i/j cols indicates that we're in this scenario
    ## - any nonzero costs will be communicated by EdgeList entries
    bookkeeping_node_labels  <- nodes %>%
        filter(is.na(upstream_not_down))$nodelabels
    bookkeeping_costs  <-
        if (any(levels(eld[['i']]) %in% bookkeeping_node_labels))
        {
                if (!anyflipped) #if `anyflipped` is T, then this 
                    eld  <- rbind(eld,#was already done previously
                                  edgelist(t(distances),
                                           node.labels(solution)
                                           )
                                  )
                bookkeeping_ij <- solution@arcs@bookkeeping %>%
                    left_join(eld, by=c("start" = "i", "end"="j")
                              )
                sum(bookkeeping_ij$dist * bookkeeping_ij$flow,
                    na.rm=TRUE)
        } else 0
    
    main_costs + bookkeeping_costs
    }
## Computing the Lagrangian given a match and a set of node prices 
##
## @param distances An InfinitySparseMatrix, DenseMatrix or EdgeList giving distances
## @param solution A MCFSolutions object
## @return The value of the Lagrangian.
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
evaluate_lagrangian <- function(distances, solution) {
    stopifnot(is(solution, "MCFSolutions"),
              is(distances, "EdgeList") ||
              is(distances, "InfinitySparseMatrix") ||
              is(distances, "DenseMatrix")
              )
    anyflipped  <- any(solution@subproblems[["flipped"]])
    ## according to Bertsekas *Network Optimization*, page 155, the Lagrangian is given by:
    ## L(x, p) = \sum_{i,j} x_{ij} (a_ij - (p_i - p_j)) + \sum_i s_i p_i
    ## where
    ##  - x_ij is the amount of flow along ij
    ##  - a_ij is the cost of the edge ij
    ##  - p_i is the cost of node ni
    ##  - s_i is the amount of flow entering or leaving the system at i
    ##

    ## note to self, need to know if problem was flipped to get distances out of the ISM.

    nodes  <- as(nodeinfo(solution), "tbl_df")
    main_ij <- left_join(solution@arcs@matches,
                         dplyr::filter(nodes, upstream_not_down),
                         by = c("groups", "upstream" = "nodelabels")) %>%
               left_join(y = dplyr::filter(nodes, !upstream_not_down),
                         by = c("groups", "downstream" = "nodelabels"),
                         suffix = c(x = ".i", y = ".j")
                         )

    eld <- edgelist(distances, node.labels(solution))
    if (anyflipped)
        eld  <- rbind(eld, edgelist(t(distances), node.labels(solution)))

    main_ij <- left_join(main_ij, eld,
                         by = c("upstream" = "i", "downstream"= "j"),
                         suffix = c(x = "", y = ".dist")
                         )

    bookkeeping_ij <- solution@arcs@bookkeeping %>%
        left_join(nodes,
                  by = c("groups", "start" = "nodelabels")
                  ) %>%
        left_join(dplyr::filter(nodes,#assumes bookkeeping arcs terminate...
                                is.na(upstream_not_down)),#...only in bookkeeping nodes
                  by = c("groups", "end" = "nodelabels"),
                  suffix = c(x = ".i", y = ".j"))

    sum_supply_price <- sum(nodes$supply * nodes$price, na.rm=TRUE)

    sum_main_flow_cost <- sum(main_ij$dist - (main_ij$price.i - main_ij$price.j))
    bookkeeping_node_labels  <- nodes %>%
        filter(is.na(upstream_not_down))$nodelabels
    if (any(levels(eld[['i']]) %in% bookkeeping_node_labels))
        warning("Distances involving bookkeeping nodes ignored/treated as 0")
    sum_bookkeeping_flow_cost  <- 
        sum(bookkeeping_ij$flow * (0 - (bookkeeping_ij$price.i - bookkeeping_ij$price.j)))

    return(sum_main_flow_cost + sum_bookkeeping_flow_cost + sum_supply_price)
}


## Compute dual functional from distance, MCF problem description
##
## Both `solution@nodes` and `solution@arcs` are used, the former
## for node prices and the latter for upper capacities of bookkeeping
## arcs.
##
## @param distances An InfinitySparseMatrix, DenseMatrix or EdgeList giving distances
## @param solution A MCFSolutions object 
## @return Value of the dual functional, a numeric of length 1.
#' @importFrom dplyr left_join
#' @importFrom dplyr inner_join 
#' @importFrom dplyr filter
evaluate_dual <- function(distances, solution) {
    stopifnot(is(solution, "MCFSolutions"),
              is(solution, "FullmatchMCFSolutions") ||
              all(nodeinfo(solution)[is.na(nodeinfo(solution)$upstream_not_down),
                                     "name"] %in% c('(_Sink_)', '(_End_)')
                  ),
              is(distances, "EdgeList") ||
              is(distances, "InfinitySparseMatrix") ||
              is(distances, "DenseMatrix")
              )
    anyflipped  <- any(solution@subproblems[["flipped"]])
    ## according to Bertsekas *Network Optimization*, page 156-7,
    ## the dual functional is given by:
    ##
    ## The dual functional, defined on pp. 155-6 of same ref., is
    ## Q(p)   = \sum_{i,j} q(a_ij, c_ij; p_i, p_j) + \sum_i s_i p_i
    ##  where
    ##  - p_i is the price/potential of node ni
    ##  - s_i is the amount of flow entering or leaving the system at i
    ##  - u_ij is the upper capacity of edge ij
    ##  - each edge ij is taken to have lower capacity 0
    ##  - a_ij is the cost of the edge ij
    ##   - q(a_ij, c_ij; p_i, p_j) =
    ##       case_when(p_i  > a_ij + p_j ~ (a_ij + p_j - p_i)*u_ij),
    ##                 p_i <= a_ij + p_j ~ 0)
    ## This Q(p) being what you get if minimize the Lagrangian over x's
    ##  respecting capacity but not conservation of flow constraints.
    ##
    nodes  <- as(nodeinfo(solution), "tbl_df")    
    sum_supply_price <- sum(nodes$supply * nodes$price, na.rm=TRUE)

    ## calculate costs from bookkeeping edges
    ##
    bookkeeping_ij <- left_join(solution@arcs@bookkeeping,
                                nodes, 
                                by = c("groups", "start" = "nodelabels")) %>%
        left_join(y = dplyr::filter(nodes,#assumes bookkeeping arcs terminate ...
                                    is.na(upstream_not_down)),#... only in bookkeeping nodes
                  by = c("groups", "end" = "nodelabels"), 
                  suffix = c(x = ".i", y = ".j"))

    nonpositive_flowcosts_bookkeeping  <-
        pmin(0,
             bookkeeping_ij$price.j - bookkeeping_ij$price.i
             ) * bookkeeping_ij$capacity

    eld <- edgelist(distances, node.labels(solution))
    bookkeeping_node_labels  <- nodes %>%
        filter(is.na(upstream_not_down))$nodelabels
    if (any(levels(eld[['i']]) %in% bookkeeping_node_labels))
        warning("Distances involving bookkeeping nodes ignored/treated as 0")

    if (anyflipped)
        eld  <- rbind(eld, edgelist(t(distances), node.labels(solution)))
    ## now check if any treated/upstream nodes are being added; if so, bail
    ## (don't currently know how to impute prices for upstream nodes.
    ##  nor do we have logic with which to impute their supplies.)
    if (any(upstream_NA  <- is.na(nodes[['price']]) &
                !is.na(nodes[['upstream_not_down']]) &
                nodes[['upstream_not_down']]
            )
        ) {
        if (any(nodes[['name']][upstream_NA] %in%
                as.character(c(eld[['i']], eld[['j']]))
                )
            )
        stop("Cannot impute node price for upstream nodes (usually treatment) that were not included in original matching problem.")
    }
    
    ## if we've gotten this far, a missing node price means that it is a down stream node
    ## and the missing price is the lesser of the sink and the overflow bookkeeping nodes
    if (any(downstream_NA  <- is.na(nodes[['price']]) &
                !is.na(nodes[['upstream_not_down']]) &
                !nodes[['upstream_not_down']]
            )
        ) {
        for (gg in levels(factor(nodes$groups)))
            {
        price_imputation <-
            min(nodes$price[nodes$groups==gg &
                            is.na(nodes$upstream_not_down)]
                )
        nodes[downstream_NA & nodes$groups==gg,
              "price"]  <- price_imputation
        }
    }

    matchable_ij <-  eld %>% 
        inner_join(y = filter(nodes, upstream_not_down),
                  by = c("i" = "nodelabels")
                  ) %>%
        inner_join(y = filter(nodes, !upstream_not_down),
                  by = c("j" = "nodelabels"),
                  suffix = c(x =".i", y =".j")
                  ) 

    nonpositive_flowcosts_matchables <-
        pmin(0,
             matchable_ij$dist -
             (matchable_ij$price.i - matchable_ij$price.j)
             )

    return(sum_supply_price +
           sum(nonpositive_flowcosts_bookkeeping) +
           sum(nonpositive_flowcosts_matchables)
           )
}
