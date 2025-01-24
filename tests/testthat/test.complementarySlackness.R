################################################################################
## Tests of functions that test for complementarySlackness
################################################################################

context("Complementary slackness, no subproblems")

make_known_optimal <- function(flipped=FALSE) {

    x <- data.frame(row.names = c("A", "B", "C", "D", "E"),
                    z = c(rep(!flipped, 2), rep(flipped, 3)),
                    y = c(0, 2, -3, 0, 6))
    m <- match_on(z ~ y, data = x, method = "euclidean")
    m  <- as.InfinitySparseMatrix(m)
    nodes <- new("NodeInfo",
                 data.frame(stringsAsFactors = FALSE,
                     name = c("A", "B", "C", "D", "E", '(_End_)'),
                     price = c(3, 4.75, 0, 2.75, 0.25, 0),
                     upstream_not_down = c(TRUE, TRUE, FALSE, FALSE, FALSE, NA),
                     supply = c(1L, 1L, 0L, 0L, 0L, -2L),
                     groups = factor(rep("a", 6))))
    node.labels(nodes)  <- nodes[['name']]

    arcs <- new("ArcInfo",
                matches = data.frame(
                    groups = factor(rep("a", 2)),
                    upstream   = factor(c("A", "B" ), levels=node.labels(nodes)),
                    downstream = factor(c("D", "E"), levels=node.labels(nodes))
                ),
                bookkeeping = data.frame(
                    groups = factor(rep("a", 3)),
                    start = factor(c("C", "D", "E"), levels=node.labels(nodes)),
                    end = factor(c("(_End_)", "(_End_)", "(_End_)"), levels=node.labels(nodes)),
                    flow = c(0L, 1L, 1L),
                    capacity = rep(1L, 3)
                ))
    subprob  <- new("SubProbInfo",
                    data.frame(groups="a", flipped=flipped, hashed_dist=character(1),
                               resolution=NA_real_, primal_value=NA_real_, dual_value=NA_real_,
                               feasible=NA, exceedance=NA_real_, solver = NA, stringsAsFactors=FALSE)
                    )
    mcf_solution  <- new("MCFSolutions", subproblems=subprob, nodes=nodes, arcs=arcs)
    list(x = x, m = m, mcf = mcf_solution)
}

make_caliper_example <- function(flipped=FALSE) {
    opt <- make_known_optimal(flipped=flipped)
    opt$cal <- caliper(opt$m, 3, values = TRUE)
    tmp <- opt$mcf@nodes
    tmp$price <- c(3.25, 5.75, 0.25, 3.25, NA_real_, 0)
    opt$mcf@nodes <- new("NodeInfo", tmp)
    node.labels(opt$mcf)  <- setNames(tmp$name, tmp$name)

    opt$mcf@arcs <- new("ArcInfo",
                    matches = data.frame(
                        groups = factor(rep("a", 2)),
                        upstream   = factor(c("A", "B" ),levels=node.labels(opt$mcf)),
                        downstream = factor(c("C", "D"),levels=node.labels(opt$mcf))
                    ),
                    bookkeeping = data.frame(
                        groups = factor(rep("a", 2)),
                        start = factor(c("C", "D"),levels=node.labels(opt$mcf)),
                        end = factor(c("(_End_)", "(_End_)"),levels=node.labels(opt$mcf)),
                        flow = c(1L, 1L),
                        capacity = rep(1L, 2)))
    return(opt)
}

make_known_optimal_fullm <- function(flipped=FALSE)
{
    x <- data.frame(row.names = c("a", "b", "c", "d", "e"),
                    z = c(rep(!flipped, 3), rep(flipped, 2)),
                    y = c(0, 2, 4, 0, 6))
    m <- match_on(z ~ y, data = x, method = "euclidean")
    m  <- as.InfinitySparseMatrix(m)
    nodes <- new("NodeInfo",
                 data.frame(stringsAsFactors = FALSE,
                     name = c("a", "b", "c", "d", "e", '(_Sink_)', '(_End_)'),
                     price = c(6, 4, 4, 0, 0, 0, 0),
                     upstream_not_down = c(TRUE, TRUE, TRUE, FALSE, FALSE, NA, NA),
                     supply = c(2L, 2L, 2L, 0L, 0L, -2L, -4L),
                     groups = factor(rep("b", 7))))
    arcs <- new("ArcInfo",
                matches=data.frame(
                    groups=factor("b"),
                    upstream=factor(c(1, 2, 3),levels=node.labels(nodes)),
                    downstream=factor(c(4, 4, 5),levels=node.labels(nodes))
                ),
                bookkeeping=data.frame(
                    groups=factor("b"),
                    start=factor(c(1, 2, 3, 4, 5, 4, 5),levels=node.labels(nodes)),
                    end=factor(c(7, 7, 7, 7, 7, 6, 6),levels=node.labels(nodes)),
                    flow=as.integer(c(1, 1, 1, 1, 0, 1, 1)),
                    capacity=as.integer(c(1, 1, 1, 2, 2, 1, 1))
                )
                )
subprob  <- new("SubProbInfo",
                    data.frame(groups='b', flipped=flipped, hashed_dist=character(1),
                               resolution=NA_real_, primal_value=NA_real_, dual_value=NA_real_,
                               feasible=NA, exceedance=NA_real_, solver = NA, stringsAsFactors=FALSE)
                    )
    mcf_solution  <- new("MCFSolutions", subproblems=subprob, nodes=nodes, arcs=arcs)
    node.labels(mcf_solution)  <- names(node.labels(mcf_solution)) # for alignment w/ m
    list(x = x, m = m, mcf = mcf_solution)
}
test_that("Compute primal", {
    opt <- make_known_optimal()
    expect_equal(evaluate_primal(opt$m, opt$mcf), 4)
  ## repeat with a dense density matrix
    expect_equal(evaluate_primal(new("DenseMatrix", as.matrix(opt$m)), opt$mcf), 4)

    ## now do it with a calipered version of the problem.
    cal <- make_caliper_example()
    expect_equal(evaluate_primal(cal$cal, cal$mcf), 5)

    ## simple f.m. problem (with both End and Sink bookkeeping nodes)
    optfm  <- make_known_optimal_fullm()
    expect_equal(evaluate_primal(optfm$m, optfm$mcf), 4)
    ## 'flipped' variant
    opt.f  <- make_known_optimal(flipped=TRUE)
    expect_equal(evaluate_primal(opt.f$m, opt.f$mcf), 4)

})
test_that("Compute Lagrangian", {
    opt <- make_known_optimal()
    ## since the above arcs represents the optimal, the lagrangian at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_lagrangian(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_lagrangian(new("DenseMatrix", as.matrix(opt$m)), opt$mcf), 4)

    ## now do it with a calipered version of the problem.
    cal <- make_caliper_example()
    expect_equal(evaluate_lagrangian(cal$cal, cal$mcf), 5)

    ## this shouldn't change when we use the full version, since none of those edges contribute flow or capacity
    expect_equal(evaluate_lagrangian(cal$m, cal$mcf), 5)

    ## simple f.m. problem (with both End and Sink bookkeeping nodes)
    optfm  <- make_known_optimal_fullm()
    expect_equal(evaluate_lagrangian(optfm$m, optfm$mcf), 4)
    ## 'flipped' variant
    opt.f  <- make_known_optimal(flipped=TRUE)
    expect_equal(evaluate_lagrangian(opt.f$m, opt.f$mcf), 4)
})


test_that("Compute dual functional", {
    opt <- make_known_optimal()

    ## since the above arcs represent the optimum, the dual functional at this point should be equal to the
    ## primal objective function (ie., the sum of matched distances).
    expect_equal(evaluate_dual(opt$m, opt$mcf), 4)

    ## repeat with a dense density matrix
    expect_equal(evaluate_dual(new("DenseMatrix", as.matrix(opt$m)), opt$mcf), 4)

    ## now do it with a calipered version of the problem.
    cal <- make_caliper_example()
    expect_equal(evaluate_dual(cal$cal, cal$mcf), 5)

    ## the caliper isn't optimal for the full problem, so the dual decreases using the full distance matrix
    expect_equal(evaluate_dual(cal$m, cal$mcf), 2.75)

    ## simple f.m. problem (with both End and Sink bookkeeping nodes)
    optfm  <- make_known_optimal_fullm()
    expect_equal(evaluate_dual(optfm$m, optfm$mcf), 4)
    ## 'flipped' variant
    opt.f  <- make_known_optimal(flipped=TRUE)
    expect_equal(evaluate_dual(opt.f$m, opt.f$mcf), 4)
})

context("Complementary slackness, multiple subproblems")

make_known_2subprobs  <- function(flipped=c(FALSE, FALSE))
    {
        o1  <- make_known_optimal(flipped=flipped[1])
        o2  <- make_known_optimal_fullm(flipped=flipped[2])
        stopifnot(length(intersect(rownames(o1$m), rownames(o2$m)))==0,
                  length(intersect(colnames(o1$m), colnames(o2$m)))==0)
        newx  <- rbind(o1$x, o2$x)
        newm  <- matrix(Inf,
                        nrow=(nrow(o1$m)+nrow(o2$m)),
                        ncol=(ncol(o1$m)+ncol(o2$m)),
                        dimnames=list(c(rownames(o1$m), rownames(o2$m)),
                                      c(colnames(o1$m), colnames(o2$m))
                                      )
                        )
        newm[rownames(o1$m), colnames(o1$m)]  <- o1$m
        newm[rownames(o2$m), colnames(o2$m)]  <- o2$m
        newm  <- as.InfinitySparseMatrix(newm)
        newmcf  <- c(o1$mcf, o2$mcf)
        list(x=newx, m=newm, mcf=newmcf)
        }

test_that("Compute primal", {
    opt <- make_known_2subprobs()
    expect_equal(evaluate_primal(opt$m, opt$mcf), 8)
    ## partly 'flipped' variant
    opt.f  <- make_known_2subprobs(flipped=c(FALSE, TRUE))
    expect_equal(evaluate_primal(opt.f$m, opt.f$mcf), 8)
})

test_that("Compute Lagrangian", {
    opt <- make_known_2subprobs()
    ## since the above arcs represents the optimal, the lagrangian
    ## at this point should be equal to the primal
    ## objective function (ie., the sum of matched distances).
    expect_equal(evaluate_lagrangian(opt$m, opt$mcf), 8)
    ## 'flipped' variant
    opt.f  <- make_known_2subprobs(flipped=c(FALSE, TRUE))
    expect_equal(evaluate_lagrangian(opt.f$m, opt.f$mcf), 8)
})

test_that("Compute dual functional", {
    opt <- make_known_2subprobs()
    ## since the above arcs represents the optimal, the dual
    ## functional at this point should be equal to the primal
    ## objective function (ie., the sum of matched distances).
    expect_equal(evaluate_dual(opt$m, opt$mcf), 8)
    ## 'flipped' variant
    opt.f  <- make_known_2subprobs(flipped=c(FALSE,TRUE))
    expect_equal(evaluate_dual(opt.f$m, opt.f$mcf), 8)
})

context("Solvers give back optimal solutions")

test_that("Verifying solvers get correct node prices", {
  mytol <- .Machine$double.eps^(1/4)
  x <- 1:5
  z <- c(1, 1, 0, 0, 1)
  units <- paste0("u", 1:5)
  names(x) <- names(z) <- units
  df <- data.frame(x, z)
  mm <- match_on(z ~ x, data = df)

  mmm <- as.matrix(mm)
  min_dist <- sum(mmm[1:2, 1]) + mm[3, 2]

  match_lemon <- fullmatch(mm, data = df, solver = "LEMON")
  mcfs_lemon <- attr(match_lemon, "MCFSolutions")

  expect_equal(evaluate_primal(mm, mcfs_lemon), min_dist)

  expect_equal(evaluate_dual(mm, mcfs_lemon), min_dist)

  if (requireNamespace("rrelaxiv", quietly = TRUE)) {
    match_relax <- fullmatch(mm, data = df, solver = "RELAX-IV", tol=mytol/10)
    mcfs_relax <- attr(match_relax, "MCFSolutions")
    expect_equal(evaluate_primal(mm, mcfs_relax), min_dist)
    expect_equal(evaluate_dual(mm, mcfs_relax), min_dist, tolerance = mytol)
  }


  ## some examples from the nuclearplants data set

  data(nuclearplants)
  npm <- match_on(pt ~ . - pt, data = nuclearplants)
  np_lemon <- fullmatch(npm, data = nuclearplants, solver = "LEMON",
                        tol = mytol/10)
  primal_lemon <- evaluate_primal(npm, attr(np_lemon, "MCFSolutions"))
  dual_lemon <- evaluate_dual(npm, attr(np_lemon, "MCFSolutions"))

  expect_equal(primal_lemon, dual_lemon, tol = mytol)

  if (requireNamespace("rrelaxiv", quietly = TRUE)) {
    np_relax <- fullmatch(npm, data = nuclearplants, solver = "RELAX-IV",
                        tol = mytol/10)
    primal_relax <- evaluate_primal(npm, attr(np_relax, "MCFSolutions"))
    dual_relax <- evaluate_dual(npm, attr(np_relax, "MCFSolutions"))

    expect_equal(primal_relax, primal_lemon, tol = mytol)
    expect_equal(primal_relax, dual_relax, tol = mytol)
  }

  npm2 <- match_on(pr ~ cost + t1, data = nuclearplants)
  np2_lemon <- fullmatch(npm2, min.controls = 1, max.controls = 3,
                         data = nuclearplants, solver = "LEMON",
                        tol = mytol/10)
  primal2_lemon <- evaluate_primal(npm2, attr(np2_lemon, "MCFSolutions"))
  dual2_lemon <- evaluate_dual(npm2, attr(np2_lemon, "MCFSolutions"))

  expect_equal(primal2_lemon, dual2_lemon, tol = mytol)


  if (requireNamespace("rrelaxiv", quietly = TRUE)) {
    np2_relax <- fullmatch(npm2, min.controls = 1, max.controls = 3,
                           data = nuclearplants, solver = "RELAX-IV",
                        tol = mytol/10)
    primal2_relax <- evaluate_primal(npm2, attr(np2_relax, "MCFSolutions"))
    dual2_relax <- evaluate_dual(npm2, attr(np2_relax, "MCFSolutions"))

    expect_equal(primal2_relax, primal2_lemon, tol = mytol)
    expect_equal(primal2_relax, dual2_relax, tol = mytol)
  }

  lemons <-  c("CycleCancelling", "CapacityScaling",
                               "CostScaling", "NetworkSimplex")

  for (solver in lemons) {
    np2_solver <- fullmatch(npm2, min.controls = 1, max.controls = 3, data = nuclearplants, solver = LEMON(solver),
                        tol = mytol/10)
    mcf2_solver <- attr(np2_solver, "MCFSolutions")
    expect_equal(evaluate_dual(npm2, mcf2_solver), dual2_lemon, label = solver, tol = mytol)
  }

  # this data frme is just to silence some warnings from fullmatch
  d <- data.frame(rep(1, 6))
  rownames(d) <- LETTERS[1:6]

  # this little example was leading to different prices right out of fmatch
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  mm <- as.InfinitySparseMatrix(m)
  mm_lemon <- fullmatch(mm, data = d,
                        tol = mytol/10)
  mm_mcfs  <- attr(mm_lemon, "MCFSolutions")
  mm_dual <- evaluate_dual(mm, mm_mcfs)

  for (solver in lemons) {
    np2_solver <- fullmatch(mm, solver = LEMON(solver), data = d,
                        tol = mytol/10)
    mcf2_solver <- attr(np2_solver, "MCFSolutions")
    expect_equal(evaluate_dual(mm, mcf2_solver), mm_dual, label = solver, tol = mytol)
  }

})
