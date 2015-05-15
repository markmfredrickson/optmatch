#############################################################################
# scores: wrapper of predict to automatically use lm()'s data as newdata
#############################################################################

library(testthat)

context("scores function")

test_that("Works like predict", {
  data(nuclearplants)
  pg <- lm(cost ~ . - pr, data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- scores(pg, newdata=nuclearplants)
  expect_equal(pred, scores1)

### error if missing newdata without being in a model (unlike predict)
  expect_error(scores(pg))
})

test_that("Works like predict with 'with' and 'attach'", {
  data(nuclearplants)
  pg <- lm(cost ~ . - pr, data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  scores1 <- with(nuclearplants, scores(pg))

  expect_equal(pred, scores1, check.attributes=FALSE)
})

test_that("Correct in a model", {
  data(nuclearplants)
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + scores(pg), data=nuclearplants)
  ps2 <- glm(pr ~ cap + date + t1 + bw + scores(pg, newdata=nuclearplants), data=nuclearplants)
  ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
  expect_equal(fitted(ps1), fitted(ps3), check.attributes=FALSE)
})

test_that("Correct in a model using 'with'", {
  data(nuclearplants)
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  options(warn=-1)
  ps1 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
  ps2 <- with(nuclearplants, glm(pr ~ cap + date + t1 + bw + scores(pg)))

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
})

test_that("Works in match_on", {
  data(nuclearplants)
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

test_that("Handles weights with no missingness", {

  data(nuclearplants)

  mod1 <- lm(pr ~ cost, weights=t1, data=nuclearplants)
  expect_true(length(scores(mod1, newdata=nuclearplants)) == 32)

  mod2 <- lm(pr ~ scores(mod1), data=nuclearplants)
  expect_true(length(mod2$fitted) == 32)

})

test_that("NA imputation", {
  data(nuclearplants)
  nuclearplants[c(3,5,12,19),1] <- NA
  nuclearplants$cap[7] <- NA

  g1 <- glm(pr ~ cost + t1 + t2, data=nuclearplants, family=binomial)
  l1 <- lm(pr ~ ne + scores(g1), data=nuclearplants)

  # check that all cases are used
  expect_equal(length(fitted(l1)), 32)

  # Filling data first should be equivalent
  np2 <- fill.NAs(pr ~ cap + cost + t1 + t2 + ne, data=nuclearplants)
  g2 <- glm(pr ~ . - ne, data=np2, family=binomial)
  l2 <- lm(pr ~ ne + scores(g2), data=np2)
  l3 <- lm(pr ~ ne + predict(g2, newdata=np2), data=np2)
  expect_equal(fitted(l2), fitted(l3))
  # But it should be different because the first stage had missingingness
  expect_false(all(fitted(l1) ==  fitted(l3)))

  # if we fill.NAs twice, should be the same.
  g1.filled <- glm(pr ~ cost + t1 + t2 + cost.NA, data=np2, family=binomial)
  l1.filled <- lm(pr ~ ne + scores(g1.filled), data=np2)
  expect_equal(fitted(l1), fitted(l1.filled))

  # Subset shouldn't affect imputation
  pgscore <- lm(cost + t1 ~ cap + interaction(ct,bw) + t2, subset=(pr==0),
                data=nuclearplants)

  l5 <- lm(pr ~ cap + scores(pgscore), data=nuclearplants)
  l6 <- lm(pr ~ cap + scores(pgscore, newdata=nuclearplants), data=nuclearplants)
  l7 <- lm(pr ~ cap + predict(pgscore, newdata=nuclearplants), data=nuclearplants)

  np3 <- fill.NAs(cost ~ cap + pr + bw + ct + t1 + t2, data=nuclearplants)
  pgscore.filled <- lm(cost + t1 ~ cap + interaction(ct,bw) + t2 + cap.NA,
                       subset=(pr==0), data=np3)
  l8 <- lm(pr ~ cap + scores(pgscore.filled), data=np3)
  l9 <- lm(pr ~ cap + scores(pgscore.filled, newdata=np3), data=np3)
  l10 <- lm(pr ~ cap + predict(pgscore.filled, newdata=np3), data=np3)

  expect_equal(fitted(l5), fitted(l6))
  expect_equal(fitted(l5), fitted(l7))
  expect_equal(fitted(l8), fitted(l9))
  expect_equal(fitted(l8), fitted(l10))
  # we don't expect l5 and l8 to be the same because np3 mean imputed over the
  # entire data set, whereas scores(pgscore) only mean imputed over the subset
  expect_equal(length(fitted(l5)), 31)
  expect_equal(length(fitted(l8)), 32)

  pgscore <- lm(cost~cap + interaction(ct,bw), subset=(pr==0),
                weights=t1, data=nuclearplants)

  w1 <- lm(pr ~ t2 + scores(pgscore, newdata=nuclearplants), data=nuclearplants)
  w2 <- lm(pr ~ t2 + scores(pgscore), data=nuclearplants)
  w3 <- lm(pr ~ t2 + predict(pgscore, newdata=nuclearplants), data=nuclearplants)
  expect_equal(length(fitted(w1)), 32)
  expect_equal(fitted(w1), fitted(w2))

  expect_equal(length(fitted(w3)), 31)
  expect_false(fitted(w1)[1] == fitted(w3)[1])
})
