#############################################################################
# Tests that should not be pushed to CRAN
#############################################################################

context("Not for CRAN")

test_that("scores with bigglm", {
  # Currently, scores() for bigglm is not supported
  if (requireNamespace("biglm", quietly = TRUE)) {
    require(biglm)
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2))

    bgps <- biglm::bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    m1 <- lm(Z ~ scores(bgps), data=test.data)
    m2 <- lm(Z ~ predict(bgps, newdata=test.data), data=test.data)

    expect_equal(fitted(m1), fitted(m2))

    test.data2 <- test.data
    test.data2$X1[1] <- NA

    # predict's na.action=na.pass is ignored in bigglm, so we'll need to impute
    # beforehand to get the same results
    bgps2 <- biglm::bigglm(Z ~ X1 + X2, data = test.data2, family = binomial())
    expect_warning(m3 <- lm(Z ~ scores(bgps2), data=test.data2),
                   "Imputation and refitting of bigglm objects")
    m4 <- lm(Z ~ predict(bgps2, newdata=fill.NAs(test.data2)), data=test.data2)

    expect_equal(fitted(m3), fitted(m4))

    suppressWarnings(bgps3 <- biglm::bigglm(Z ~ X1*X2, data = test.data2, family = binomial()))
    fill.test.data <- fill.NAs(Z ~ X1*X2, data=test.data2)
    expect_warning(m5 <- lm(Z ~ scores(bgps3), data=fill.test.data),
                   "Imputation and refitting of bigglm objects")
    m6 <- lm(Z ~ predict(bgps3, newdata=fill.test.data), data=test.data2)

    expect_equal(fitted(m5), fitted(m6))
  }
  expect_true(TRUE) # avoid empty test warning
})

test_that("match_on with bigglm distances", {
  if (requireNamespace("biglm", quietly = TRUE)) {
    require(biglm)
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))


    bgps <- biglm::bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    res.bg <- match_on(bgps, data = test.data)

    # compare to glm
    res.glm <- match_on(glm(Z ~ X1 + X2, data = test.data, family = binomial()))
    expect_equivalent(res.bg, res.glm)
  }
  expect_true(TRUE) # avoid empty test warning
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


  if (requireNamespace("rrelaxiv", quietly = TRUE)) {
  # Not necessarily seeing the same speedups on LEMON, so only running this for
  # RELAX-IV.
#  for (i in 1:2) {
#    if (i == 1 & requireNamespaceN"amespace", quietly = TRUE("rrelaxiv", quietly = TRUE)) {
    slvr <- "RELAX-IV"
#    } else {
#      slvr <- "LEMON"
#    }

    t0  <- system.time(res0 <- fmatch(pm, 1, 1, 1, f=ff,
                                      node_info=nodes_dummy, solver = slvr)
                       )
    expect_false(is.null(mcfs0  <-  res0$MCFSolution))
    n0  <-  mcfs0@nodes
    t1  <- system.time(res1 <- fmatch(pm, 1, 1, 1, f=ff, node_info=n0,
                                      solver = slvr))
    expect_gt(t0['elapsed'],
              t1['elapsed'])
    #  }
  }
  expect_true(TRUE) # avoid empty test warning
})
### next tests disabled due to an odd scoping-related error
### occurring only within the test, not in interactive use.
### (The problem is that `test.design` can't be found when
### the model is re-fit inside of `scores()`.)
###test_that("scores with svyglm (survey package) objects",
if (FALSE) {
  if (requireNamespace("survey", quietly = TRUE)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))
    test.design <- svydesign(~1, probs=1, data=test.data)
    sglm <- svyglm(Z ~ X1 + X2, test.design, family = quasibinomial())
    expect_silent(scores(sglm, newdata=sglm$data))
}
}
###)
###test_that("match_on with svyglm (survey package) objects",
if (FALSE) {
  if (requireNamespace("survey", quietly = TRUE)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2),
                            B = rep(c(0,1), times = n/2))
    test.design <- svydesign(~1, probs=1, data=test.data)
    sglm <- svyglm(Z ~ X1 + X2, test.design, family = binomial())
    expect_silent(res.svy0 <- match_on(sglm, data=test.data, standardization.scale=1))
    expect_silent(res.svy1 <- match_on(sglm, data=test.data, standardization.scale=svy_sd))
    expect_silent(res.svy2 <- match_on(sglm, data=test.data, standardization.scale=svy_mad))
    ## Comparisons to glm: currently failing, disabled pending investigation.
    ## First step: figure out whether sglm/aglm are returning the same
    ## linear predictors.
    if (FALSE)
    {
        aglm <- glm(Z ~ X1 + X2, test.data, family = binomial())
    res.glm0 <- match_on(aglm, data=test.data, standardization.scale=1)
        res.glm1 <- match_on(aglm, data=test.data, standardization.scale=stats::sd)
        mad_lo  <- function(x) {
            ctr  <- quantile(x, probs=0.5, type=1)
            stats::mad(x, center=ctr, low=TRUE)
                       }
    res.glm2 <- match_on(aglm, data=test.data, standardization.scale=mad_lo)

    expect_equivalent(res.svy0, res.glm0)
    expect_equivalent(res.svy1, res.glm1)
    expect_equivalent(res.svy2, res.glm2)
    }
}
}
###)


test_that("survival::strata masking doesn't break", {
  data(nuclearplants)

  # survey depends on survival, so need to remove both prior to this test
  try(detach("package:survey"), silent = TRUE)
  try(detach("package:survival"), silent = TRUE)

  expect_true(!any(grepl("package:survival", search())))

  m1 <- match_on(pr ~ cost, within=exactMatch(pr ~ pt, data=nuclearplants),
                  data=nuclearplants)
  m2 <- match_on(pr ~ cost + strata(pt), data=nuclearplants)
  m2b <- match_on(pr ~ cost, data=nuclearplants)

  expect_true(is(m1, "BlockedInfinitySparseMatrix"))

  expect_true(is(m1, "BlockedInfinitySparseMatrix"))

  expect_true(is(m2, "BlockedInfinitySparseMatrix"))
  expect_true(all.equal(m1, m2, check.attributes=FALSE))
  expect_true(!isTRUE(all.equal(m2, m2b, check.attributes=FALSE)))

  if (requireNamespace("survival", quietly = TRUE)) {
    require(survival) # requireNamespace doesn't load into search path, so we
                      # need to actually load it here

    # Since we have our own optmatch::survival, a basic test
    # to ensure things still work if **survival** gets loaded and
    # survival::strata overloads optmatch::strata

    expect_true(any(grepl("package:survival", search())))

    sm2 <- match_on(pr ~ cost + strata(pt), data=nuclearplants)

    expect_identical(m2, sm2)

    detach("package:survival")

  }



})


test_that("cox model testing", {
  # this is the test case from:
  # https://github.com/markmfredrickson/optmatch/issues/44

  if (requireNamespace("survival", quietly = TRUE)) {
    data(heart, package = "survival")
    coxps <- predict(survival::coxph(survival::Surv(start, stop, event) ~
                                       age + year + transplant + cluster(id),
                                     data=heart))
    names(coxps) <- row.names(heart)
    coxmoA <- match_on(coxps, z = heart$event, caliper = 1)
    expect_true(max(coxmoA) <= 1)

    coxmoC <- match_on(coxps, within = exactMatch(event ~ transplant,
                                                  data = heart),
                       z = heart$event, caliper = 1)
    expect_true(max(coxmoC) <= 1)
  }
})
