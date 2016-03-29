
context("fullmatch function old")

test_that("fullmatch function, old version", {
  data(plantdist)
  # this will give a warning about not having data order
  expect_warning(f <- fullmatch(1 * (plantdist < 10)), "Without 'data' argument the order of the match is not guaranteed
    to be the same as your original data.") # make plantdist < 10 numeric, not logical
  expect_equal(as.vector(summary.factor(f)), c(2,2,2,3,2,2,13))
  expect_equal(names(summary.factor(f)), c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7"))

  data(nuclearplants)
  mhd2 <- match_on(pr ~ date + cum.n, data = nuclearplants,
                   within = exactMatch(pr ~ pt, data = nuclearplants))
  # the previous version of optmatch used fullmatch(mhd2 < 1)
  # this is the equivalent using an ISM (logical operators treat them as numeric
  # vectors)
  mhd2[mhd2 < 1] <- 1
  mhd2[mhd2 >= 1] <- 0
  g <- fullmatch(mhd2, data = nuclearplants)
  expect_equal(as.vector(summary.factor(g)), c(2,2,2,2,2,2,14,2,2,2))
  expect_equal(names(summary.factor(g)), c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "1.1", "1.2", "1.3"))

})
