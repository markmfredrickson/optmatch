#############################################################################
# scores: wrapper of predict to automatically use lm()'s data as newdata
#############################################################################

library(testthat)

context("scores function")

test_that("Works like predict", {
  options(warn=-1)
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- scores(pg, newdata=nuclearplants)
  expect_equal(pred, scores1)

### error if missing newdata without being in a model (unlike predict)
  expect_error(scores(pg))
})

test_that("Works like predict with 'with'", {
  options(warn=-1)
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))
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
