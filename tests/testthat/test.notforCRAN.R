#############################################################################
# Tests that should not be pushed to CRAN
#############################################################################

context("Not for CRAN")

test_that("scores with bigglm", {
  # Currently, scores() for bigglm is not supported
  if (require(biglm)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2))

    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    m1 <- lm(Z ~ scores(bgps), data=test.data)
    m2 <- lm(Z ~ predict(bgps, newdata=test.data), data=test.data)

    expect_equal(fitted(m1), fitted(m2))

    test.data2 <- test.data
    test.data2$X1[1] <- NA

    # predict's na.action=na.pass is ignored in bigglm, so we'll need to impute
    # beforehand to get the same results
    bgps2 <- bigglm(Z ~ X1 + X2, data = test.data2, family = binomial())
    expect_warning(m3 <- lm(Z ~ scores(bgps2), data=test.data2),
                   "Imputation and refitting of bigglm objects")
    m4 <- lm(Z ~ predict(bgps2, newdata=fill.NAs(test.data2)), data=test.data2)

    expect_equal(fitted(m3), fitted(m4))

    suppressWarnings(bgps3 <- bigglm(Z ~ X1*X2, data = test.data2, family = binomial()))
    fill.test.data <- fill.NAs(Z ~ X1*X2, data=test.data2)
    expect_warning(m5 <- lm(Z ~ scores(bgps3), data=fill.test.data),
                   "Imputation and refitting of bigglm objects")
    m6 <- lm(Z ~ predict(bgps3, newdata=fill.test.data), data=test.data2)

    expect_equal(fitted(m5), fitted(m6))
  }

})

test_that("match_on with bigglm distances", {
  if (require(biglm)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))


    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    res.bg <- match_on(bgps, data = test.data)

    # compare to glm
    res.glm <- match_on(glm(Z ~ X1 + X2, data = test.data, family = binomial()))
    expect_equivalent(res.bg, res.glm)
  }
})

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
test_that("Hinting decreases runtimes",{
  v <- c(1, Inf, 2,
         2, 1, Inf,
         3, 2, 1)
  suppressWarnings(v  <- as.integer(v))
  # the clear match to make: 
  # A:D, B:E, C:F
  m <- matrix(v, nrow = 3, ncol = 3)
  colnames(m) <- c("A", "B", "C")
  rownames(m) <- c("D", "E", "F")
  pm <- edgelist(m)

  t0  <- system.time(res <- fmatch(pm, 2, 2,
                                   node_info=pairmatch_nodeinfo(pm))
                     )
  expect_false(is.null(mcfs0  <-  res$MCFSolution))
  n0  <-  mcfs0@nodes
  t1  <- system.time(fmatch(pm, 2, 2, node_info=n0))
  expect_gt(t0['elapsed'], #to be honest I don't consistently see a difference here. To avoid failure
            t1['elapsed']) #we might set this to an `expect_gte`, although it would be better to fix it.
})
