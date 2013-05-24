#############################################################################
# scores: wrapper of predict to auto-matically use lm()'s data as newdata
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

test_that("Works like predict with 'with' and 'attach'", {
  options(warn=-1)
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- with(nuclearplants, scores(pg))
  attach(nuclearplants)
  scores2 <- scores(pg)
  detach(nuclearplants)
  expect_equal(pred, scores1, check.attributes=FALSE)
  expect_equal(pred, scores2, check.attributes=FALSE)

})

test_that("Correct in a model", {
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + scores(pg), data=nuclearplants)
  ps2 <- glm(pr ~ cap + date + t1 + bw + scores(pg, newdata=nuclearplants), data=nuclearplants)
  ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)

# Check all comparisons (ignore 22:24, which will differ due to formula naming)
  expect_equal(ps1[-c(22:24)], ps2[-c(22:24)], check.attributes=FALSE)
  expect_equal(ps1[-c(22:24)], ps3[-c(22:24)], check.attributes=FALSE)
})

test_that("Correct in a model using 'with' or 'attach'", {
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
  ps2 <- with(nuclearplants, glm(pr ~ cap + date + t1 + bw + scores(pg)))
  attach(nuclearplants)
  expect_error(ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg)))
  # Am I doing something wrong in how attach would be used with predict here?
  # Note that while predict bombs in this scenario, scores actually works!
  ps4<- glm(pr ~ cap + date + t1 + bw + scores(pg))
  detach(nuclearplants)

  expect_equal(ps1[-c(22:25)], ps2[-c(22:25)], check.attributes=FALSE)
  expect_equal(ps1[-c(22:25)], ps4[-c(22:25)], check.attributes=FALSE)
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

