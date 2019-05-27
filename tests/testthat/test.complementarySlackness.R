################################################################################
## Tests of functions that test for complementarySlackness
################################################################################

context("Complementary slackness")

test_that("Compute Lagrangian", {
    ##
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

    ## since the above arcs represents the optimal, the lagrangian at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_lagrangian(m, nodes, arcs), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_lagrangian(as.matrix(m), nodes, arcs), 4)
})


test_that("Compute dual functional", {
    ##
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
 
    ## since the above arcs represent the optimum, the dual functional at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_dual(m, nodes, arcs), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_dual(as.matrix(m), nodes, arcs), 4)
})
