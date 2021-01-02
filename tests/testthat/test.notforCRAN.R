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
  N <- 1000
  X <- data.frame(X1 = rnorm(N), 
                  X2 = rnorm(N, mean = runif(N, -5, 5)), 
                  X3 = as.factor(sample(letters[1:5], N, replace = T)))
  mm <- model.matrix(I(rep(1, N)) ~  X1 + X2 + X1:X3, data = X)
  coefs <- runif(dim(mm)[2], -2, 2)
  logits <- as.vector(coefs %*% t(mm)) 
  DATA <- data.frame(Z = rbinom(N, size = 1, prob = plogis(logits)), X)  
  m <- match_on(x = Z ~ X1 + X2 + X3, data = DATA)
  if (nrow(m) > ncol(m)) m <- t(m)
  ff <- nrow(m)/ncol(m)
  pm <- edgelist(m)
  pm$dist <- as.integer(100*pm$dist)
  nodes_dummy <- pairmatch_nodeinfo(pm)
  t0  <- system.time(res0 <- fmatch(pm, 1, 1, 1, f=ff,
                                   node_info=nodes_dummy)
                     )
  expect_false(is.null(mcfs0  <-  res0$MCFSolution))
  n0  <-  mcfs0@nodes
  t1  <- system.time(res1 <- fmatch(pm, 1, 1, 1, f=ff, node_info=n0))
  expect_gt(t0['elapsed'], 
            t1['elapsed'])
})
