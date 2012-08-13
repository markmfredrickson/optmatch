################################################################################
# Tests for interaction with RItools
################################################################################

if(require("testthat") && require("RItools")) {
context("RItools interaction")

test_that("Summary function adds RItools info", {

  # fill.NAs is a layer of indirection between the glm object and RItools
  # putting th test here because this could happen for other ways to create glms, not just fill.NAs
  test.data <- data.frame(Z = rep(c(0,1),10), X = c(1,2,NA,3,4,NA,5,NA,9,10))
  test.glm <- glm(fill.NAs(Z ~ X, data = test.data), family = binomial)

  test.pm <- pairmatch(mdist(test.glm))

  # may cause an error if the function can't get at the proper data that was used to create the glm object
  res <- summary(test.pm, test.glm)

  expect_true(!is.null(res$balance))
  
})

}
