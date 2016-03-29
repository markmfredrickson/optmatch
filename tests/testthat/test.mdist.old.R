
context("mdist function old")

test_that("Basic test", {
  data(nuclearplants)
  fmla <- pr ~ t1 + t2 + pt
  test.glm <- glm(fmla, family=binomial(), data=nuclearplants)

  result.glm <- mdist(test.glm)
  result.glm2 <- mdist(test.glm, ~pt)

  expect_true(inherits(result.glm, "optmatch.dlist"))
  expect_true(length(result.glm) == 1)
  expect_true(length(result.glm2) == 2)

  result.optmatch.dlist <- mdist(result.glm)
  expect_true(inherits(result.optmatch.dlist, "optmatch.dlist"))
  expect_true(length(result.optmatch.dlist) == 1)

  result.fmla <- mdist(fmla, data = nuclearplants)
  expect_true(inherits(result.fmla, "optmatch.dlist"), "Should be a optmatch object")
  expect_true(length(result.fmla) == 1)

  result.fmla2 <- mdist(fmla, ~pt, data = nuclearplants)
  expect_true(inherits(result.fmla2, "optmatch.dlist"), "Should be a optmatch object")
  expect_true(length(result.fmla2) == 2)

  ### Can the results be updated?
  # update for optmatch.dlist objects is failing
  #expect_true(identical(result.glm2, update(result.glm, structure.fmla=~pt)))
  #expect_true(identical(result.fmla2, update(result.fmla, structure.fmla=~pt)))
})

test_that("Function Tests", {
  data(nuclearplants)
  fmla <- pr ~ t1 + t2 + pt
  result.fmla <- mdist(fmla, data = nuclearplants)
  result.fmla2 <- mdist(fmla, ~pt, data = nuclearplants)
  test.glm <- glm(fmla, family=binomial(), data=nuclearplants)

  # first, a simpler version of scalar diffs
  sdiffs <- function(treatments, controls) {
    abs(outer(treatments$t1, controls$t1, `-`))
  }

  result.function <- mdist(sdiffs, pr ~ 1, nuclearplants)
  expect_equal(dim(result.function), c(10,22))

  expect_equal(optmatch:::parseFmla(y ~ a | group), lapply(c("y", "a", "group"), as.name))
  expect_equal(optmatch:::parseFmla(y ~ a), c(as.name("y"), as.name("a"), NULL))

  expect_true(!is.null(rownames(result.function$m)) && all(rownames(result.function$m) %in% rownames(nuclearplants[nuclearplants$pr == 1,])))

  stripCall <- function(obj) {
    attr(obj, "call") <- NULL
    obj
  }

  result.function.a <- mdist(sdiffs, pr ~ 1 | pt, nuclearplants)
  result.function.b <- mdist(sdiffs, pr ~ pt, nuclearplants)

  expect_equal(stripCall(result.function.a), stripCall(result.function.b))

  ### Check of updating:
#  expect_equal(stripCall(result.function.b), stripCall(update(result.function,structure.fmla=pr~pt)))

  expect_error(mdist(sdiffs, pr ~ pt + t1, nuclearplants))

  # the fun part, making a dlist when there are multiple groups

  expect_equal(length(mdist(sdiffs, pr ~ pt, nuclearplants)), 2)

  ### Using mad() instead of sd() for GLM distances

  result <- mdist(glm(pr ~ t1 + t2 + cost, data = nuclearplants, family = binomial()))

  # this is an odd test, but a simple way to make sure mad is running, not SD().
  # I would like a better test of the actual values, but it works
  expect_true(mean(result$m) > 2)

  ### mdist() should informatively complain if passed a numeric vector.
  ### (This may change in the future.)

  expect_error(mdist(test.glm$linear.predictor))

  ### Stratifying by a pipe (|) character in formulas

  main.fmla <- pr ~ t1 + t2
  strat.fmla <- ~ pt
  combined.fmla <- pr ~ t1 + t2 | pt

  result.main <- mdist(main.fmla, structure.fmla = strat.fmla, data = nuclearplants)
  result.combined <- mdist(combined.fmla, data = nuclearplants)

  expect_equal(stripCall(result.main), stripCall(result.combined))

  ### Informatively insist that one of formulas specify the treatment group
  expect_error(mdist(~t1+t2, structure.fmla=~pt, data=nuclearplants))
  expect_equal(stripCall(mdist(pr~t1+t2, structure.fmla=~pt, data=nuclearplants)),
                 stripCall(mdist(~t1+t2, structure.fmla=pr~pt, data=nuclearplants)))

  ### Finding "data" when it isn't given as an argument
  ### Caveats:
  ### * data's row.names get lost when you don't pass data as explicit argument;
  ### thus testing with 'all.equal(unlist(<...>),unlist(<...>))' rather than 'identical(<...>,<...>)'.
  ### * with(nuclearplants, mdist(fmla)) bombs for obscure scoping-related reasons,
  ### namely that the environment of fmla is the globalenv rather than that created by 'with'.
  ### This despite the facts that identical(fmla,pr ~ t1 + t2 + pt) is TRUE and that
  ### with(nuclearplants, mdist(pr ~ t1 + t2 + pt)) runs fine.
  ### But then with(nuclearplants, lm(fmla)) bombs too, for same reason, so don't worry be happy.
  attach(nuclearplants)
  expect_equal(unlist(result.fmla),unlist(mdist(fmla)))
  expect_equal(unlist(result.main),unlist(mdist(main.fmla, structure.fmla=strat.fmla)))
  expect_equal(unlist(result.combined),unlist(mdist(combined.fmla)))
  detach("nuclearplants")
  expect_equal(fmla,pr ~ t1 + t2 + pt)
  expect_equal(unlist(result.fmla),unlist(with(nuclearplants, mdist(pr ~ t1 + t2 + pt))))
  expect_equal(combined.fmla, pr ~ t1 + t2 | pt)
  expect_equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2 | pt))))
  expect_equal(unlist(result.fmla), unlist(with(nuclearplants[-which(names(nuclearplants)=="pt")],
                                                  mdist(update(pr ~ t1 + t2 + pt,.~.-pt + nuclearplants$pt))
                                                  )
                                             )
                 )
  expect_equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2, structure.fmla=strat.fmla))))
})

test_that('bigglm method', {
  data(nuclearplants)

  fmla <- pr ~ t1 + t2 + pt
  result.fmla <- mdist(fmla, data = nuclearplants)
  result.fmla2 <- mdist(fmla, ~pt, data = nuclearplants)

  test.glm <- glm(fmla, family=binomial(), data=nuclearplants)
  result.glm <- mdist(test.glm)
  result.glm2 <- mdist(test.glm, ~pt)

  # first, a simpler version of scalar diffs
  sdiffs <- function(treatments, controls) {
    abs(outer(treatments$t1, controls$t1, `-`))
  }

  result.function <- mdist(sdiffs, pr ~ 1, nuclearplants)
  result.function.a <- mdist(sdiffs, pr ~ 1 | pt, nuclearplants)

  if (require('biglm')) {
    bgps <- bigglm(fmla, data=nuclearplants, family=binomial() )
    expect_error(mdist(bgps, structure.fmla=pr ~ 1))
    expect_error(mdist(bgps, data=nuclearplants))
    result.bigglm1 <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants)
    result.bigglm1a <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
                             standardization.scale=sd)
    result.bigglm1b <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
                             standardization.scale=NULL)
    result.bigglm2 <- mdist(bgps, structure.fmla=pr ~ pt, data=nuclearplants)
  }

  ### Jake found a bug 2010-06-14
  ### Issue appears to be a missing row.names/class

  absdist1 <- mdist(sdiffs, structure.fmla = pr ~ 1|pt, data = nuclearplants)
  expect_warning(pabsdist1 <- pairmatch(absdist1))
  expect_true(length(pabsdist1) > 0)

  ### Check that distances combine as they should
  ### (a joint test of mdist and Ops.optmatch.dlist)
  ### Distances without subclasses:
  expect_true(inherits(result.glm + result.fmla, "optmatch.dlist"))
  expect_true(inherits(result.glm + result.function, "optmatch.dlist"))
  if (require("biglm"))
    expect_true(inherits(result.glm + result.bigglm1, "optmatch.dlist"))


  ### Distances embodying subclassification:
  expect_true(inherits(result.glm2 + result.fmla2, "optmatch.dlist"))

  expect_true(inherits(result.glm2 + result.function.a, "optmatch.dlist"))


  if (require("biglm"))
    expect_true(inherits(result.glm2 + result.bigglm2, "optmatch.dlist"))
})
