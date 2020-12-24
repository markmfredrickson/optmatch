################################################################################
### R/Fortran Interface Tests
################################################################################

context("R/Fortran Interface")

# convenience function for use in testing
pairmatch_nodeinfo  <- function(edges) {
    stopifnot(is(edges, "EdgeList"))
    allunits  <- levels(edges[['i']])
    istreated  <- allunits %in% edges[['i']]

    adf  <- data.frame(name=c(allunits, "(_Sink_)"),
                       price=0L,
                       upstream_not_down=c(istreated, NA),
                       supply=c(rep(1L, sum(istreated)),
                                rep(0L, sum(!istreated)),
                                -sum(istreated)
                                ),
                       stringsAsFactors=FALSE
                       )
    new("NodeInfo", adf)
}

test_that("fmatch accepts DistanceSpecifications", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  # and doesn't accept other things...
  expect_error(fmatch(v, 2, 2))

  # the goal of this matrix is that there is a clear match to make
  # A:D, B:E, C:F
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  pm <- edgelist(m)

  res <- fmatch(pm, 2, 2, node_info=pairmatch_nodeinfo(pm))
  expect_true(all(c("j","i",
                    "dist", # used in `doubleSolve()`'s "maxerr" calc
                    "solution") %in% names(res)))
  expect_equal(length(res$solution), 7) # seven non-Inf entries

  # check that A-D is a pair and A-B is not a match
  expect_equal(res$solution[res$j == "A" & res$i == "D"], 1)
  expect_equal(res$solution[res$j == "A" & res$i == "B"],
               numeric(0))

  M <- as.InfinitySparseMatrix(m)
  pM <- edgelist(M)
  res.ism <- fmatch(pM, 2, 2, node_info=pairmatch_nodeinfo(pM))
  expect_identical(res$solution, res.ism$solution)
})

test_that("Stop on unacceptable input", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")  
  rownames(m) <- c("D", "E", "F")  
  
  m1  <- m
  colnames(m1) <- c("(_Sink_)", "B", "C")  
  pm1  <- edgelist(m1)
  expect_error(fmatch(pm1,2,2, node_info=pairmatch_nodeinfo(pm1)), "unique") #"(_Sink_)"

  m2  <- m1
  colnames(m2) <- c("A", "B", "C")
  rownames(m2) <- c("(_End_)", "E", "F")

  pm2  <- edgelist(m2)
  expect_error(fmatch(pm2,2,2, node_info=pairmatch_nodeinfo(pm2)), "(_End_)")

})

test_that("Solutions -> factor helper", {
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)

  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")

  skeleton <- edgelist(m)
  class(skeleton)  <- "data.frame" #drops S4 class
  skeleton  <- dplyr::mutate(skeleton, treated=factor(i), control=factor(j))

  pairs <- cbind(skeleton, solution = c(1,0,0,1,0,0,1))
  pairs.expected <- c(1,2,3,1,2,3)
  names(pairs.expected) <- c("D", "E", "F", "A", "B", "C")

  expect_equal(solution2factor(pairs), pairs.expected)

  pairOfTriples <- cbind(skeleton, solution = c(1,0,1,0,0,1,1))
  pot.expected <- c(1,2,2,1,1,2)
  names(pot.expected) <- c("D", "E", "F", "A", "B", "C")
  expect_equal(solution2factor(pairOfTriples), pot.expected)

  treatedNotMatched <- cbind(skeleton, solution = c(1,0,0,1,1,0,0))
  tnm.expected <- c(1,2, NA, 1,2,1)
  names(tnm.expected) <- c("D", "E", "F", "A", "B", "C")

  expect_equal(solution2factor(treatedNotMatched), tnm.expected)

  controlNotMatched <- cbind(skeleton, solution = c(0,0,1,1,0,0,1))
  cnm.expected <- c(1, 1, 3, NA, 1, 3)
  names(cnm.expected) <- c("D", "E", "F", "A", "B", "C")

  expect_equal(solution2factor(controlNotMatched), cnm.expected)

  # handles failed matchings by returning NULL
  noMatches <- cbind(skeleton, solution = -1)

  expect_true(is.null(solution2factor(noMatches)))
})  

test_that("Passing and receiving node information",{
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)
  # the clear match to make: 
  # A:D, B:E, C:F
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  pm <- edgelist(m)

  res <- fmatch(pm, 2, 2, node_info=pairmatch_nodeinfo(pm))
  expect_false(is.null(mcfs0  <-  res$MCFSolution))
  n0  <-  mcfs0@nodes
  expect_silent(fmatch(pm, 2, 2, node_info=n0))
  n0_madebad  <- n0
  expect_is(n0_madebad$price, "integer")
  n0_madebad[n0_madebad$name=="A", 'price']  <- .5 # no longer integer
  expect_error(fmatch(pm, 2, 2, node_info=n0_madebad))

  expect_false(n0[n0$name=="A",'upstream_not_down']) # 'A' is downstream,
  n1  <- new("NodeInfo", n0[n0$name!="A",])#  so we can pass a
  expect_gt(nrow(n0), nrow(n1)) # NodeInfo that doesn't mention it.
  expect_silent(fmatch(pm, 2, 2, node_info=n1))
})
