################################################################################
# Fullmatch-recover-infeasible tests
################################################################################

context("fullmatch-recover-infeasible update")

# basic tests are in the general test.fullmatch.R file (as this update should not change
# functionality for fundamentally feasible problems)

test_that("Invalid mean.controls input", {
  data(nuclearplants)
  m <- match_on(glm(pr~cost, data=nuclearplants))
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = 1, omit.fraction=.5),
               "omit.fraction and mean.controls cannot both be specified")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = -1),
               "mean.controls must be NULL or numeric greater than 0")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = "a"),
               "mean.controls must be NULL or numeric greater than 0")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = 23/10),
               "mean.controls cannot be larger than the ratio of number of controls to treatments")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = c(1,1)))
               # "Length of 'mean.controls' arg must be same as number of subproblems [1]")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = 1, min.controls=2),
               "mean.controls cannot be smaller than min.controls")
  expect_error(fullmatch(m, data=nuclearplants, mean.controls = 1, max.controls=1/2),
               "mean.controls cannot be larger than max.controls")

  set.seed(2)
  x <- runif(20)
  fact <- c(rep(0,7), rep(1, 4), rep(2, 9))
  treat <- c(rep(0,4), rep(1, 2),0, rep(0, 2), rep(1, 2), rep(0, 5), rep(1, 4))
  dd <- as.data.frame(cbind(x,fact,treat))

  mm <- match_on(treat~.-fact, data=dd, within=exactMatch(treat~fact, dd))

  expect_error(fullmatch(mm, data=dd, mean.controls=1.5),
               "mean.controls cannot be larger than the ratio of number of controls to treatments")
  expect_error(fullmatch(mm, data=dd, mean.controls=c(2, NA, 2)),
               "mean.controls cannot be larger than the ratio of number of controls to treatments")

})

test_that("mean.controls should do the same as omit.fraction", {
  data(nuclearplants)
  m <- match_on(glm(pr~cost, data=nuclearplants))
  f <- fullmatch(m, data=nuclearplants, omit.fraction=1/2)
  g <- fullmatch(m, data=nuclearplants, mean.controls=11/10)

  expect_true(all.equal(f, g, check.attributes=FALSE))

  set.seed(2)
  x <- runif(20)
  fact <- c(rep(0,7), rep(1, 4), rep(2, 9))
  treat <- c(rep(0,4), rep(1, 2),0, rep(0, 2), rep(1, 2), rep(0, 5), rep(1, 4))
  dd <- as.data.frame(cbind(x,fact,treat))

  mm <- match_on(treat~.-fact, data=dd, within=exactMatch(treat~fact, dd))

  f <- fullmatch(mm,data=dd, omit.fraction=c(1/5, 1/2, 1/5))
  g <- fullmatch(mm,data=dd, mean.controls=c(2, 1/2, 1))

  expect_true(all.equal(f, g, check.attributes=FALSE))

  f <- fullmatch(mm,data=dd, omit.fraction=c(1/5, NA, NA))
  g <- fullmatch(mm,data=dd, mean.controls=c(2, NA, NA))

  expect_true(all.equal(f, g, check.attributes=FALSE))

  f <- fullmatch(mm,data=dd, omit.fraction=c(3/5, NA, 1/5))
  g <- fullmatch(mm,data=dd, mean.controls=1)

  expect_true(all.equal(f, g, check.attributes=FALSE))

})

test_that("Allow passing all NA to mean.controls or omit.fraction", {
  data(nuclearplants)
  m <- exactMatch(pr ~ pt, data=nuclearplants)

  f <- fullmatch(m, data=nuclearplants)

  fm1 <- fullmatch(m, data=nuclearplants, mean.controls=c(NA,NA))
  fm2 <- fullmatch(m, data=nuclearplants, mean.controls=NA)
  fm3 <- fullmatch(m, data=nuclearplants, mean.controls=c(NULL, NULL))
  fm4 <- fullmatch(m, data=nuclearplants, mean.controls=NULL)
  fm5 <- fullmatch(m, data=nuclearplants, mean.controls=c(NA, NULL))

  fo1 <- fullmatch(m, data=nuclearplants, omit.fraction=c(NA,NA))
  fo2 <- fullmatch(m, data=nuclearplants, omit.fraction=NA)
  fo3 <- fullmatch(m, data=nuclearplants, omit.fraction=c(NULL, NULL))
  fo4 <- fullmatch(m, data=nuclearplants, omit.fraction=NULL)
  fo5 <- fullmatch(m, data=nuclearplants, omit.fraction=c(NA, NULL))

  attr(f, "call") <- NULL
  attr(fm1, "call") <- attr(fm2, "call") <- attr(fm3, "call") <- attr(fm4, "call") <- attr(fm5, "call") <- NULL
  attr(fo1, "call") <- attr(fo2, "call") <- attr(fo3, "call") <- attr(fo4, "call") <- attr(fo5, "call") <- NULL

  expect_true(identical(f, fm1))
  expect_true(identical(f, fm2))
  expect_true(identical(f, fm3))
  expect_true(identical(f, fm4))
  expect_true(identical(f, fm5))

  expect_true(identical(f, fo1))
  expect_true(identical(f, fo2))
  expect_true(identical(f, fo3))
  expect_true(identical(f, fo4))
  expect_true(identical(f, fo5))

})

test_that("Correctly apply max.controls", {
  options("optmatch_verbose_messaging" = TRUE)

  set.seed(2)
  x <- runif(20)
  fact <- c(rep(0,7), rep(1, 4), rep(2, 9))
  treat <- c(rep(0,4), rep(1, 2),0, rep(0, 2), rep(1, 2), rep(0, 5), rep(1, 4))
  dd <- as.data.frame(cbind(x,fact,treat))

  mm <- match_on(treat~.-fact, data=dd, within=exactMatch(treat~fact, dd))

  # have three subgroups:
  # 1) 5 ctrl, 2 treat
  # 2) 2 ctrl, 2 treat
  # 3) 5 ctrl, 4 treat

  # no restrictions, everything should be matched.
  s1 <- stratumStructure(f <- fullmatch(mm,data=dd))
  expect_true(all(unlist(strsplit(names(s1), ":")) > 0))

  # max.controls = 2
  expect_warning(s2 <- stratumStructure(g <- fullmatch(mm,data=dd, max.controls=2)))
  max.controls <- max(as.numeric(unlist(lapply(strsplit(names(s2), ":"),"[",2))))
  expect_true(max.controls <= 2)

  # max controls = 1
  expect_warning(s3 <- stratumStructure(h <- fullmatch(mm,data=dd, max.controls=1)))
  max.controls <- max(as.numeric(unlist(lapply(strsplit(names(s3), ":"),"[",2))))
  expect_true(max.controls <= 1)

  # size of control group is sum of treatment group of pmin of
  # max.controls and control:treatment ratio for tx group member
  # (prior to resolution of issue 74, the below led to a single 2:1
  # matched set)
  adist <- matrix(c(1:4, rep(Inf, 8)), 2, 6, dimnames=list(letters[1:2], letters[3:8]))
  expect_silent(fullmatch(adist, data=data.frame(1:8, row.names=letters[1:8])))
  expect_warning(fm <- fullmatch(adist, max.c=1, data=data.frame(1:8, row.names=letters[1:8])), "infeasible")

  expect_true(all(table(fm)==2))
})

test_that("Omits occur only on controls", {
  options("optmatch_verbose_messaging" = TRUE)
  set.seed(3)
  x <- runif(20)
  fact <- c(rep(0,7), rep(1, 4), rep(2, 9))
  treat <- c(rep(0,4), rep(1, 2),0, rep(0, 2), rep(1, 2), rep(0, 5), rep(1, 4))
  dd <- as.data.frame(cbind(x,fact,treat))

  mm <- match_on(treat~.-fact, data=dd, within=exactMatch(treat~fact, dd))

  expect_warning(s1 <- stratumStructure(fullmatch(mm, data=dd, max.controls=2)))
  ctrls1 <- as.numeric(unlist(lapply(strsplit(names(s1), ":"),"[",2)))
  treats1 <- as.numeric(unlist(lapply(strsplit(names(s1), ":"),"[",1)))
  # It should drop some of the controls (some treats1 should be 0)
  # but none of the treatments (all ctrls1 > 0)
  expect_true(all(ctrls1 > 0))
  expect_true(any(treats1 == 0))

  expect_warning(s2 <- stratumStructure(fullmatch(mm, data=dd, max.controls=1)))
  ctrls2 <- as.numeric(unlist(lapply(strsplit(names(s2), ":"),"[",2)))
  treats2 <- as.numeric(unlist(lapply(strsplit(names(s2), ":"),"[",1)))
  expect_true(all(ctrls2 > 0))
  expect_true(any(treats2 == 0))
})

test_that("If omit.fraction is included", {
  set.seed(10)
  x <- runif(20)
  fact <- c(rep(0,7), rep(1, 4), rep(2, 9))
  treat <- c(rep(0,4), rep(1, 2),0, rep(0, 3), rep(1, 1), rep(0, 5), rep(1, 4))
  dd <- as.data.frame(cbind(x,fact,treat))

  mm <- match_on(treat~.-fact, data=dd, within=exactMatch(treat~fact, dd))

  # have three subgroups:
  # 1) 5 ctrl, 2 treat
  # 2) 3 ctrl, 1 treat
  # 3) 5 ctrl, 4 treat

  f <- fullmatch(mm,data=dd, omit.fraction=c(1/5, 1/3, 1/5))
  # check that exactly 1 is omitted from each.
  expect_true(sum(is.na(f[row.names(dd[dd$fact == 0 & dd$treat == 0,])])) == 1)
  expect_true(sum(is.na(f[row.names(dd[dd$fact == 1 & dd$treat == 0,])])) == 1)
  expect_true(sum(is.na(f[row.names(dd[dd$fact == 2 & dd$treat == 0,])])) == 1)

  expect_warning(g <- fullmatch(mm,data=dd, max.controls=1, omit.fraction=c(1/5, 1/3, 1/5)))
  # infeasible even though some omit.fraction is given, needs to drop more
  expect_true(sum(is.na(g[row.names(dd[dd$fact == 0 & dd$treat == 0,])])) == 3)
  expect_true(sum(is.na(g[row.names(dd[dd$fact == 1 & dd$treat == 0,])])) == 2)
  expect_true(sum(is.na(g[row.names(dd[dd$fact == 2 & dd$treat == 0,])])) == 1)
})

test_that("Suggested omit.fraction can be used", {
  data(nuclearplants)

  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants)

  expect_warning(s1 <- fullmatch(mm, data=nuclearplants, max.controls=2))
  s2 <- fullmatch(mm, data=nuclearplants, max.controls=2, omit.fraction=attr(s1, "omit.fraction"))

  expect_warning(s3 <- fullmatch(mm, data=nuclearplants, max.controls=1))
  s4 <- fullmatch(mm, data=nuclearplants, max.controls=1, omit.fraction=attr(s3, "omit.fraction"))

  expect_true(all.equal(s1, s2, check.attributes=FALSE))
  expect_true(all.equal(s3, s4, check.attributes=FALSE))
})

test_that("mean.controls as fraction", {
  data(nuclearplants)

  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants)

  # 22 treatments, I want to exclude 3, so 19ctrls/10treat = 1.9 mean.controls
  s1 <- stratumStructure(fullmatch(mm, data=nuclearplants, mean.controls=1.9, max.controls=3))
  expect_equal(3, sum(s1[substr(names(s1),1,1) == 0]))

  # 22, exclude 10, 12/10 = 1.2 mean.controls
  s2 <- stratumStructure(fullmatch(mm, data=nuclearplants, mean.controls=1.2, max.controls=3))
  expect_equal(10, sum(s2[substr(names(s2),1,1) == 0]))
})

test_that("attr saved after recovery", {
  data(nuclearplants)

  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants)

  # not infeasible as given
  f <- fullmatch(mm, data=nuclearplants, max.controls=3)
  expect_equal(attr(f, "min.controls"), 0)
  expect_equal(attr(f, "mean.controls"), NULL)
  expect_equal(attr(f, "max.controls"), 3)
  expect_equal(attr(f, "omit.fraction"), as.numeric(NA))

  # infeasible as given
  expect_warning(f <- fullmatch(mm, data=nuclearplants, max.controls=2))
  s <- stratumStructure(fullmatch(mm, data=nuclearplants))
  expect_equal(attr(f, "min.controls"), 0)
  expect_equal(attr(f, "mean.controls"), NULL)
  expect_equal(attr(f, "max.controls"), 2)
  # how many SHOULD we omit?
  numomit <- sum(pmax(0,as.numeric(unlist(lapply(strsplit(names(s), "1:"), "[", 2)))-2))
  expect_equal(attr(f, "omit.fraction"), numomit/22)


  # not infeasible as mean.controls provided
  f <- fullmatch(mm, data=nuclearplants, max.controls=2, mean.controls=2)
  expect_equal(attr(f, "min.controls"), 0)
  expect_equal(attr(f, "mean.controls"), 2)
  expect_equal(attr(f, "max.controls"), 2)
  expect_equal(attr(f, "omit.fraction"), NULL)


  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants, within=exactMatch(pr ~ pt, data=nuclearplants))

  # infeasible as given for subproblem 1, feasible for subproblem 2
  expect_warning(f <- fullmatch(mm, data=nuclearplants, max.controls=2))
  expect_equal(attr(f, "min.controls"), c(0,0), check.attributes=FALSE)
  expect_equal(attr(f, "mean.controls"), NULL)
  expect_equal(attr(f, "max.controls"), c(2,2), check.attributes=FALSE)
  expect_equal(attr(f, "omit.fraction"), c(9/19, NA), check.attributes=FALSE)

  Z <- c(1,0,0,0,0,1,0,0)
  B <- c(rep('a', 5), rep('b', 3))
  d <- data.frame(Z, B)

  res.b <- exactMatch(Z ~ B, data=d)

  expect_warning(f <- fullmatch(res.b, data=d, max.controls=2))
  a <- c(0,0)
  names(a) <- c('a','b')
  expect_equal(attr(f, "min.controls"), a)
  expect_equal(attr(f, "mean.controls"), NULL)
  a <- c(2,2)
  names(a) <- c('a','b')
  expect_equal(attr(f, "max.controls"), a)
  a <- c(1/2,NA)
  names(a) <- c('a','b')
  expect_equal(attr(f, "omit.fraction"), a)


})

test_that("fullmatch_try_recovery", {
  data(nuclearplants)

  mm <- match_on(pr ~ cost + t1 + t2, data=nuclearplants)

  options("fullmatch_try_recovery" = TRUE)
  # warn and fix
  expect_warning(expect_true(any(is.na(fullmatch(mm, data=nuclearplants, max.controls = 2)))))
  options("fullmatch_try_recovery" = FALSE)
  # fail to fix
  expect_warning(expect_true(all(is.na(fullmatch(mm, data=nuclearplants, max.controls = 2)))))
  options("fullmatch_try_recovery" = TRUE)
  # back to fixing.
  expect_warning(expect_true(any(is.na(fullmatch(mm, data=nuclearplants, max.controls = 2)))))

})



test_that("n_t > n_c", {
  data(nuclearplants)

  nuclearplants$pr <- abs(1-nuclearplants$pr)
  # 22 treatment, 10 control
  m <- match_on(pr ~ cost, data=nuclearplants)

  # should pass here without problems
  expect_true(any(!is.na(fullmatch(m, data=nuclearplants))))

  # min.controls = 1/2, so we need 11 controls. Can't accomodate.
  expect_warning(expect_true(all(is.na(fullmatch(m, min.controls = 1/2, data=nuclearplants)))))
})

test_that("Issue #92", {
  # Based upon data that had 1058 controls, 62 treated, and min=1, max=5, omit=.8.
  d <- data.frame("z" <- c(rep(0,1058), rep(1, 62)),
                  "x" <- rnorm(1120))

  expect_that(fullmatch(z ~ x, data=d, min=1, max=5), gives_warning("infeasible"))
  # this shouldn't warn about infeasible
  expect_silent(fullmatch(z ~ x, data=d, min=1, max=5, omit=.8))
  expect_that(fullmatch(z ~ x, data=d, min=1, max=2, omit=.2), gives_warning("infeasible"))
})
