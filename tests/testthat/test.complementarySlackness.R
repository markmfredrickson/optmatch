################################################################################
## Tests of functions that test for complementarySlackness
################################################################################

context("Complementary slackness")

make_known_optimal <- function() {

    x <- data.frame(row.names = c("A", "B", "C", "D", "E"), z = c(1, 1, 0, 0, 0), y = c(0, 2, -3, 0, 6))
    m <- match_on(z ~ y, data = x, method = "euclidean")
    nodes <- new("NodeInfo",
                 data.frame(stringsAsFactors = FALSE,
                     name = c("A", "B", "C", "D", "E", '(_End_)'),
                     price = c(3, 4.75, 0, 2.75, 0.25, 0),
                     upstream_not_down = c(TRUE, TRUE, FALSE, FALSE, FALSE, NA),
                     supply = c(1L, 1L, 0L, 0L, 0L, -2L),
                     groups = factor(rep("a", 6))))

    arcs <- new("ArcInfo",
                matches = data.frame(
                    groups = rep("a", 2),
                    upstream   = factor(c("A", "B" )),
                    downstream = factor(c("D", "E"))),
                bookkeeping = data.frame(
                    groups = rep("a", 3),
                    start = factor(c("C", "D", "E")),
                    end = factor(c("(_End_)", "(_End_)", "(_End_)")),
                    flow = c(0L, 1L, 1L),
                    capacity = rep(1L, 3)
                ))
    subprob  <- new("SubProbInfo",
                    data.frame(groups=character(1), flipped=FALSE, hashed_dist=character(1),
                               resolution=NA_real_, lagrangian_value=NA_real_, dual_value=NA_real_,
                               feasible=NA, exceedance=NA_real_, stringsAsFactors=FALSE)
                    )
    mcf_solution  <- new("MCFSolutions", subproblems=subprob, nodes=nodes, arcs=arcs,
                         matchables=new("MatchablesInfo"))

    list(x = x, m = m, mcf = mcf_solution)
}

test_that("Compute Lagrangian", {
    opt <- make_known_optimal()
    ## since the above arcs represents the optimal, the lagrangian at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_lagrangian(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_lagrangian(as.matrix(opt$m), opt$mcf), 4)
})


test_that("Compute dual functional", {
    opt <- make_known_optimal()
 
    ## since the above arcs represent the optimum, the dual functional at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_dual(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_dual(as.matrix(opt$m), opt$mcf), 4)
})
