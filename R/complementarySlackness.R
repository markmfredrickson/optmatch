## Computing the Lagrangian given a mathc and a set of node prices 
##
## @param distances An InfinitySparseMatrix giving distances
## @param solution A MCFSolutions object
## @return The value of the Lagrangian.
#' @importFrom dplyr left_join
evaluate_lagrangian <- function(distances, solution) {
    stopifnot(is(solution, "MCFSolutions"),
              nrow(solution@subproblems)==1)
    flipped  <- solution@subproblems[1L, "flipped"]
    ## according to Bertsekas *Network Optimization*, page 155, the Lagrangian is given by:
    ## L(x, p) = \sum_{i,j} x_{ij} (a_ij - (p_i - p_j)) + \sum_i s_i p_i
    ## where
    ##  - x_ij is the amount of flow along ij
    ##  - a_ij is the cost of the edge ij
    ##  - p_i is the cost of node ni
    ##  - s_i is the amount of flow entering or leaving the system at i
    ##

    ## note to self, need to know if problem was flipped to get distances out of the ISM.

    main_ij <- left_join(solution@arcs@matches,
                         subset(solution@nodes, upstream_not_down),
                         by = c("upstream" = "name")) %>%
               left_join(y = subset(solution@nodes, !upstream_not_down),
                         by = c("downstream" = "name"),
                         suffix = c(x = ".i", y = ".j"))

    eld <- edgelist(distances)
    if (!flipped) {
        main_ij <- left_join(main_ij,
                             eld,
                             by = c("upstream" = "i", "downstream"= "j"),
                             suffix = c(x = "", y = ".dist"))
    } else {
        main_ij <- left_join(main_ij,
                             eld,
                             by = c("upstream" = "j", "downstream"= "i"),
                             suffix = c(x = "", y = ".dist"))
    }

    bookkeeping_ij <- left_join(solution@arcs@bookkeeping,
                                as.data.frame(unclass(solution@nodes)),
                                by = c("start" = "name")) %>%
        left_join(y = subset(solution@nodes,
                             is.na(upstream_not_down)),#assumes bookkeeping arcs...
                  by = c("end" = "name"),#...terminate only in bookkeeping nodes
                  suffix = c(x = ".i", y = ".j"))


    sum_supply_price <- sum(solution@nodes$supply * solution@nodes$price)

    sum_flow_cost <- sum(main_ij$dist - (main_ij$price.i - main_ij$price.j)) +
        sum(bookkeeping_ij$flow * (0 - (bookkeeping_ij$price.i - bookkeeping_ij$price.j)))

    return(sum_flow_cost + sum_supply_price)
}


## Compute dual functional from distance, MCF problem description
##
## Both `solution@nodes` and `solution@arcs` are used, the former
## for node prices and the latter for upper capacities of bookkeeping
## arcs.
##
## @param distances An InfinitySparseMatrix giving distances
## @param solution A MCFSolutions object 
## @return Value of the dual functional, a numeric of length 1.
#' @importFrom dplyr left_join
evaluate_dual <- function(distances, solution) {
    stopifnot(is(solution, "MCFSolutions"),
              nrow(solution@subproblems)==1)
    flipped  <- solution@subproblems[1L, "flipped"]
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
    sum_supply_price <- sum(solution@nodes$supply * solution@nodes$price)

    ## calculate costs from bookkeeping edges
    ##
    bookkeeping_ij <- left_join(solution@arcs@bookkeeping,
                                as.data.frame(unclass(solution@nodes)),
                                by = c("start" = "name")) %>%
        left_join(y = subset(solution@nodes,
                             is.na(upstream_not_down)), #assumes bookkeeping arcs...
                             by = c("end" = "name"), #... terminate only in bookkeeping nodes
                             suffix = c(x = ".i", y = ".j"))

    nonpositive_flowcosts_bookkeeping  <-
        pmin(0,
             bookkeeping_ij$price.j - bookkeeping_ij$price.i
             ) * bookkeeping_ij$capacity

    ## now do edges corresponding to potential matches
    ##
    eld <- edgelist(distances)
    ## this time we have to pay attn to whether problem was flipped
    ##
    suffices  <-
        if (!flipped) c(x =".i", y =".j") else c(x =".j", y =".i")
    matchable_ij <- left_join(eld,
                         subset(solution@nodes, upstream_not_down),
                         by = c("i" = "name")) %>%
               left_join(y = subset(solution@nodes, !upstream_not_down),
                         by = c("j" = "name"),
                         suffix = suffices) # if nec., flip right at end.

    nonpositive_flowcosts_matchables <-
        pmin(0,
             matchable_ij$dist - (matchable_ij$price.i - matchable_ij$price.j)
             )

    return(sum_supply_price +
           sum(nonpositive_flowcosts_bookkeeping) +
           sum(nonpositive_flowcosts_matchables))
}
