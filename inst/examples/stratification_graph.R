data(nuclearplants)
### Full matching on a Mahalanobis distance.
fm1 <- fullmatch(pr ~ t1 + t2, data = nuclearplants)
g0 <- stratification_graph(fm1, node_color = nuclearplants$pr) + ggtitle("Full MH Stratification")
## Uncomment this next to see the graph
## print(g0)
