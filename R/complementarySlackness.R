## Computing the Lagrangian given a mathc and a set of node prices 
##
## @param distances An InfinitySparseMatrix giving distances
## @param nodes A NodeInfo object (a data.frame with specific columns)
## @param arcs A ArcInfo object (two data frames that hold the matches and the booking arcs)
## @return The value of the Lagrangian.
evaluate_lagrangian <- function(distances, nodes, arcs, flipped = FALSE) {
    ## according to Bertsekas *Network Optimization*, page 155, the Lagrangian is given by:
    ## L(x, p) = \sum_{i,j} x_{ij} (a_ij - (p_i - p_j)) + \sum_i s_i p_i
    ## where
    ##  - x_ij is the amount of flow along ij
    ##  - a_ij is the cost of the edge ij
    ##  - p_i is the cost of node ni
    ##  - s_i is the amount of flow entering or leaving the system at i
    ##

    ## note to self, need to know if problem was flipped to get distances out of the ISM.

    main_ij <- left_join(arcs@matches,
                         nodes,
                         by = c("upstream" = "name")) %>%
               left_join(y = nodes,
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

    bookkeeping_ij <- left_join(arcs@bookkeeping,
                                nodes,
                                by = c("start" = "name")) %>%
                      left_join(y = nodes,
                                by = c("end" = "name"),
                                suffix = c(x = ".i", y = ".j"))


    sum_supply_price <- sum(nodes$supply * nodes$price)

    sum_flow_cost <- sum(main_ij$dist - (main_ij$price.i - main_ij$price.j)) +
        sum(bookkeeping_ij$flow * (0 - (bookkeeping_ij$price.i - bookkeeping_ij$price.j)))

    return(sum_flow_cost + sum_supply_price)
}
