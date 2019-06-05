################################################################################
## Tests of functions that test for complementarySlackness
################################################################################

context("Complementary slackness")

make_known_optimal <- function() {

    x <- data.frame(row.names = c("A", "B", "C", "D", "E"), z = c(1, 1, 0, 0, 0), y = c(0, 2, -3, 0, 6))
    m <- match_on(z ~ y, data = x, method = "euclidean")
    nodes <- new("NodeInfo",
                 data.frame(stringsAsFactors = FALSE,
                     name = c("A", "B", "C", "D", "E", '(_End_)', '(_Sink_)'),
                     price = c(3, 4.75, 0, 2.75, 0.25, 0, 0),
                     upstream_not_down = c(TRUE, TRUE, FALSE, FALSE, FALSE, NA, NA),
                     supply = c(1L, 1L, 0L, 0L, 0L, -2L, 0L),
                     groups = factor(rep("a", 7))))

    arcs <- new("ArcInfo",
                matches = data.frame(
                    groups = rep("a", 2),
                    upstream   = factor(c("A", "B" )),
                    downstream = factor(c("D", "E"))),
                bookkeeping = data.frame(
                    groups = rep("a", 4),
                    start = factor(c("C", "D", "E", "C")),
                    end = factor(c("(_End_)", "(_End_)", "(_End_)", "(_Sink_)")),
                    flow = c(0L, 1L, 1L,0L),
                    capacity = rep(1L, 4)
                ))
    subprob  <- new("SubProbInfo",
                    data.frame(groups=character(1), flipped=FALSE, hashed_dist=character(1),
                               resolution=NA_real_, lagrangian_value=NA_real_, dual_value=NA_real_,
                               feasible=NA, exceedance=NA_real_, stringsAsFactors=FALSE)
                    )
    mcf_solution  <- new("FullmatchMCFSolutions", subproblems=subprob, nodes=nodes, arcs=arcs,
                         matchables=new("MatchablesInfo"))

    list(x = x, m = m, mcf = mcf_solution)
}

make_caliper_example <- function() {
    opt <- make_known_optimal()
    opt$cal <- caliper(opt$m, 3, values = TRUE)
    tmp <- opt$mcf@nodes[c(1:4, 6), ]
    tmp$price <- c(3.25, 5.75, 0.25, 3.25, 0)
    opt$mcf@nodes <- new("NodeInfo", tmp)

    opt$mcf@arcs <- new("ArcInfo",
                    matches = data.frame(
                        groups = rep("a", 2),
                        upstream   = factor(c("A", "B" )),
                        downstream = factor(c("C", "D"))),
                    bookkeeping = data.frame(
                        groups = rep("a", 2),
                        start = factor(c("C", "D")),
                        end = factor(c("(_End_)", "(_End_)")),
                        flow = c(1L, 1L),
                        capacity = rep(1L, 2)))
    return(opt)
}

test_that("Compute Lagrangian", {
    opt <- make_known_optimal()
    ## since the above arcs represents the optimal, the lagrangian at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_lagrangian(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_lagrangian(as.matrix(opt$m), opt$mcf), 4)

    ## now do it with a calipered version of the problem.
    cal <- make_caliper_example()
    expect_equal(evaluate_lagrangian(cal$cal, cal$mcf), 5)

    ## this shouldn't change when we use the full version, since none of those edges contribute flow or capacity
    expect_equal(evaluate_lagrangian(cal$m, cal$mcf), 5)
})


test_that("Compute dual functional", {
    opt <- make_known_optimal()
 
    ## since the above arcs represent the optimum, the dual functional at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_dual(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_dual(as.matrix(opt$m), opt$mcf), 4)

    ## now do it with a calipered version of the problem.
    cal <- make_caliper_example()
    expect_equal(evaluate_dual(cal$cal, cal$mcf), 5)

    ## the caliper isn't optimal for the full problem, so the dual decreases using the full distance matrix
    expect_equal(evaluate_dual(cal$m, cal$mcf), 2.75)
})
