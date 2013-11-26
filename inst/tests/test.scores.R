#############################################################################
# scores: wrapper of predict to automatically use lm()'s data as newdata
#############################################################################

library(testthat)

context("scores function")

test_that("Works like predict", {
  pg <- lm(cost ~ . - pr, data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- scores(pg, newdata=nuclearplants)
  expect_equal(pred, scores1)

### error if missing newdata without being in a model (unlike predict)
  expect_error(scores(pg))
})

test_that("Works like predict with 'with' and 'attach'", {
  pg <- lm(cost ~ . - pr, data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- with(nuclearplants, scores(pg))

  expect_equal(pred, scores1, check.attributes=FALSE)
})

test_that("Correct in a model", {
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + scores(pg), data=nuclearplants)
  ps2 <- glm(pr ~ cap + date + t1 + bw + scores(pg, newdata=nuclearplants), data=nuclearplants)
  ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
  expect_equal(fitted(ps1), fitted(ps3), check.attributes=FALSE)
})

test_that("Correct in a model using 'with'", {
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
  ps2 <- with(nuclearplants, glm(pr ~ cap + date + t1 + bw + scores(pg)))

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
})

test_that("Works in match_on", {
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  m1 <- match_on(pr~cap + predict(pg, newdata=nuclearplants), data=nuclearplants)
  m2 <- match_on(pr~cap + scores(pg), data=nuclearplants)
  m3 <- with(match_on(pr~cap + scores(pg)), data=nuclearplants)

  expect_equal(m1@.Data, m2@.Data)
  expect_true(all.equal(m1@.Data, m3@.Data, check.attributes=FALSE))

})

test_that("Finds variables", {
  data(iris)
  iris$Species <- as.numeric(iris$Species=="setosa")
  mod <- lm(Sepal.Length ~ . - Species, data=iris, subset=(Species==1))
  psa <- glm(Petal.Width ~ scores(mod), data=iris) # to build its own newdata, scores() is forced
  #to find a variable (Petal.Length) not OW mentioned in containing formula
  psb <- glm(Petal.Width ~ predict(mod, newdata=iris), data=iris)
  expect_equal(fitted(psa), fitted(psb))
})

test_that("indep vars of class logical properly handled", {
  data(iris)
  iris$Species <- iris$Species=="setosa"
  mod <- lm(Sepal.Length ~ Species + Sepal.Width, data=iris, subset=Species)
  psb <- glm(Petal.Width ~ Species + predict(mod, newdata=iris), data=iris)
  psa <- glm(Petal.Width ~ Species + scores(mod), data=iris)
  expect_equal(fitted(psa), fitted(psb))
})

test_that("NA imputation", {
  data(nuclearplants)
  nuclearplants[c(3,5,12,19),1] <- NA

  g <- glm(pr ~ cost + t1 + t2, data=nuclearplants, family=binomial)
  l1 <- lm(pr ~ cap + scores(g), data=nuclearplants)

  # check that all cases are used
  expect_equal(length(fitted(l1)), 32)

  # should be the same result if we fill.NA'd first.
  np2 <- fill.NAs(pr ~ cap + cost + t1 + t2, data=nuclearplants)
  g2 <- glm(pr ~ . - cap, data=np2, family=binomial)
  l2 <- lm(pr ~ cap + scores(g2), data=np2)
  l3 <- lm(pr ~ cap + predict(g2, newdata=np2), data=np2)

  # not including cost.NA should be different
  g3 <- glm(pr ~ cost + t1 + t2, data=np2, family=binomial)
  l4 <- lm(pr ~ cap + scores(g3), data=np2)

  expect_equal(fitted(l1), fitted(l2))
  expect_equal(fitted(l1), fitted(l3))
  expect_true(!all(fitted(l1) == fitted(l4)))

  pgscore <- lm(cost~cap + interaction(ct,bw), weights=t1, subset=(pr==0), data=nuclearplants)

  # This doesn't work without the newdata argument.
  expect_error(l5 <- lm(pr ~ cap + scores(pgscore), data=nuclearplants))
  l6 <- lm(pr ~ cap + scores(pgscore, newdata=nuclearplants), data=nuclearplants)
  l7 <- lm(pr ~ cap + predict(pgscore, newdata=nuclearplants), data=nuclearplants)

  np3 <- fill.NAs(pr ~ cap + cost + bw + ct, data=nuclearplants)
  l8 <- lm(pr ~ cap + scores(pgscore), data=np3)
  l9 <- lm(pr ~ cap + scores(pgscore, newdata=np3), data=np3)
  l10 <- lm(pr ~ cap + predict(pgscore, newdata=np3), data=np3)

  expect_true(!exists("l5"))
  expect_true(all.equal(fitted(l7), fitted(l6)))
  expect_true(all.equal(fitted(l7), fitted(l8)))
  expect_true(all.equal(fitted(l7), fitted(l9)))
  expect_true(all.equal(fitted(l7), fitted(l10)))
})
