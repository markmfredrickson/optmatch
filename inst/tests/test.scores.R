#############################################################################
# scores: wrapper of predict to auto-matically use lm()'s data as newdata
#############################################################################

library(testthat)

context("scores function")

test_that("Basic Tests", {
  options(warn=-1)
  data(iris)
  iris$Species <- as.numeric(as.numeric(iris$Spec) > 1)
  mod <- lm(Sepal.Length ~ ., data=iris, subset=(Species==0))
### function as predict correctly
  expect_equal(scores(mod, newdata=iris), predict(mod, newdata=iris))

### error if missing newdata (unlike predict)
  expect_error(scores(mod))
})

test_that("Correct Scores()", {
  data(iris)
  iris$Species <- as.numeric(as.numeric(iris$Spec) > 1)
  mod <- lm(Sepal.Length ~ ., data=iris, subset=(Species==0))

  options(warn=-1)
  ps1 <- glm(Petal.Width ~ Petal.Length + scores(mod), data=iris)
  ps2 <- glm(Petal.Width ~ Petal.Length + scores(mod, newdata=iris), data=iris)
  ps3 <- glm(Petal.Width ~ Petal.Length + predict(mod, newdata=iris), data=iris)

# Check all comparisons (ignore 22:24, which will differ due to formula naming)
  expect_equal(ps1[-c(22:24)], ps2[-c(22:24)], check.attributes=FALSE)
  expect_equal(ps1[-c(22:24)], ps3[-c(22:24)], check.attributes=FALSE)
  expect_equal(ps2[-c(22:24)], ps3[-c(22:24)], check.attributes=FALSE)
})

test_that("Finds variables", {
 data(iris)
 iris$Species <- as.numeric(iris$Species=="setosa")
 mod <- lm(Sepal.Length ~ ., data=iris, subset=Species)
 psa <- glm(Petal.Width ~ scores(mod), data=iris) # to build its own newdata, scores() is forced
                                        #to find a variable (Petal.Length) not OW mentioned in containing formula
 psb <- glm(Petal.Width ~ predict(mod, newdata=iris), data=iris)
 expect_equal(fitted(psa), fitted(psb))
})

### Fails. Let's fix!
###test_that("indep vars of class logical properly handled", {
### data(iris)
### iris$Species <- iris$Species=="setosa"
### mod <- lm(Sepal.Length ~ Species + Sepal.Width, data=iris, subset=Species)
### psb <- glm(Petal.Width ~ Species + predict(mod, newdata=iris), data=iris)
### psa <- glm(Petal.Width ~ Species + scores(mod), data=iris) # bombs (as of commit 060f033690d827d5944bda8836c8a76d74ac216c)
### expect_equal(fitted(psa), fitted(psb))
###})

