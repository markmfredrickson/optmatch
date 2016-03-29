#############################################################################
# scores: wrapper of predict to automatically use lm()'s data as newdata
#############################################################################

context("scores function")

test_that("Works like predict", {
  data(nuclearplants)
  pg1 <- lm(cost ~ . - pr, data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred1 <- predict(pg1, newdata=nuclearplants)
  scores1 <- scores(pg1, newdata=nuclearplants)
  expect_equal(pred1, scores1)

### error if missing newdata without being in a model (unlike predict)
  expect_error(scores(pg1))

  pg2 <- lm(pr ~ cap*cost, data=nuclearplants)
  pred2 <- predict(pg2, newdata=nuclearplants)
  scores2 <- scores(pg2, newdata=nuclearplants)

  expect_equal(pred2, scores2)

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

  suppressWarnings({
    ps1 <- glm(pr ~ cap + date + t1 + bw + scores(pg), data=nuclearplants)
    ps2 <- glm(pr ~ cap + date + t1 + bw + scores(pg, newdata=nuclearplants), data=nuclearplants)
    ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
  })

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
  expect_equal(fitted(ps1), fitted(ps3), check.attributes=FALSE)

})

test_that("Correct in a model using 'with'", {
  data(nuclearplants)
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  suppressWarnings({
    ps1 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
    ps2 <- with(nuclearplants, glm(pr ~ cap + date + t1 + bw + scores(pg)))
  })

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
})

test_that("Works in match_on", {
  data(nuclearplants)
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  expect_warning({
    m1 <- match_on(pr~cap + predict(pg, newdata=nuclearplants), data=nuclearplants)
    m2 <- match_on(pr~cap + scores(pg), data=nuclearplants)
    m3 <- with(match_on(pr~cap + scores(pg)), data=nuclearplants)
  },
  "rank-deficient"
  )

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
  suppressWarnings({
    psb <- glm(Petal.Width ~ Species + predict(mod, newdata=iris), data=iris)
    psa <- glm(Petal.Width ~ Species + scores(mod), data=iris)
  })
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
  # First, fill data using the formula, plus the extra column needed in the model
  np2 <- fill.NAs(update(formula(g1), ~ . + ne), data=nuclearplants)
  # Then refit, including the missingness indicator
  g2 <- glm(pr ~ cost + t1 + t2+ cost.NA, data=np2, family=binomial)
  l2 <- lm(pr ~ ne + scores(g2), data=np2)
  l3 <- lm(pr ~ ne + predict(g2, newdata=np2), data=nuclearplants)
  expect_equal(fitted(l2), fitted(l3))
  expect_equal(fitted(l1), fitted(l3))

  g3 <- glm(pr ~ cost*cap + t1, data=nuclearplants, family=binomial)
  np3 <- fill.NAs(pr ~ cost*cap + t1 + ne, data=nuclearplants)
  # R'll be upset if we try to include `cost:cap` as a variable
  # in a formula, so lets just fix that up.
  names(np3)[6] <- "costcap"
  g4 <- glm(pr ~ cost + cap + costcap + t1 + cost.NA + cap.NA,
            data=np3, family=binomial)

  suppressWarnings(l4 <- lm(pr ~ ne + scores(g3), data=nuclearplants))
  l4a <- lm(pr ~ ne + scores(g4), data=np3)
  l4b <- lm(pr ~ ne + predict(g4, data=np3), data=nuclearplants)

  expect_equal(fitted(l4), fitted(l4a))
  expect_equal(fitted(l4), fitted(l4b))

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

  nuclearplants$ct[1] <- NA
  pgscore <- lm(cost~cap + interaction(ct,bw), subset=(pr==0),
                weights=t1, data=nuclearplants)

  suppressWarnings({
    w1 <- lm(pr ~ t2 + scores(pgscore, newdata=nuclearplants), data=nuclearplants)
    w2 <- lm(pr ~ t2 + scores(pgscore), data=nuclearplants)
  })
  w3 <- lm(pr ~ t2 + predict(pgscore, newdata=nuclearplants), data=nuclearplants)
  expect_equal(length(fitted(w1)), 32)
  expect_equal(fitted(w1), fitted(w2))

  expect_equal(length(fitted(w3)), 30)
  expect_false(fitted(w1)[1] == fitted(w3)[1])


  pgscore <- lm(cost~cap + ct*bw, subset=(pr==0),
                weights=t1, data=nuclearplants)

  suppressWarnings({
    w4 <- lm(pr ~ t2 + scores(pgscore, newdata=nuclearplants), data=nuclearplants)
    w5 <- lm(pr ~ t2 + scores(pgscore), data=nuclearplants)
  })
  w6 <- lm(pr ~ t2 + predict(pgscore, newdata=nuclearplants), data=nuclearplants)
  expect_equal(length(fitted(w4)), 32)
  expect_equal(fitted(w4), fitted(w5))

  expect_equal(length(fitted(w6)), 30)
  expect_false(fitted(w4)[1] == fitted(w6)[1])

})

test_that("scores with bigglm", {
  # Currently, scores() for bigglm is not supported
  if (require(biglm)) {
    n <- 16
    test.data <- data.frame(Z = rep(0:1, each = n/2),
                            X1 = rnorm(n, mean = 5),
                            X2 = rnorm(n, mean = -2, sd = 2))

    bgps <- bigglm(Z ~ X1 + X2, data = test.data, family = binomial())
    m1 <- lm(Z ~ scores(bgps), data=test.data)
    m2 <- lm(Z ~ predict(bgps, newdata=test.data), data=test.data)

    expect_equal(fitted(m1), fitted(m2))

    test.data2 <- test.data
    test.data2$X1[1] <- NA

    # predict's na.action=na.pass is ignored in bigglm, so we'll need to impute
    # beforehand to get the same results
    bgps2 <- bigglm(Z ~ X1 + X2, data = test.data2, family = binomial())
    expect_warning(m3 <- lm(Z ~ scores(bgps2), data=test.data2),
                   "Imputation and refitting of bigglm objects")
    m4 <- lm(Z ~ predict(bgps2, newdata=fill.NAs(test.data2)), data=test.data2)

    expect_equal(fitted(m3), fitted(m4))

    suppressWarnings(bgps3 <- bigglm(Z ~ X1*X2, data = test.data2, family = binomial()))
    fill.test.data <- fill.NAs(Z ~ X1*X2, data=test.data2)
    expect_warning(m5 <- lm(Z ~ scores(bgps3), data=fill.test.data),
                   "Imputation and refitting of bigglm objects")
    m6 <- lm(Z ~ predict(bgps3, newdata=fill.test.data), data=test.data2)

    expect_equal(fitted(m5), fitted(m6))
  }

})

test_that("scores works with *, interaction and strata", {
  data(nuclearplants)
  g1 <- glm(pr ~ cost + interaction(ne,ct), data=nuclearplants, family=binomial)
  g2 <- glm(pr ~ cost + strata(ne,ct), data=nuclearplants, family=binomial)
  g3 <- glm(pr ~ cost + ne*ct, data=nuclearplants, family=binomial)

  expect_equal(predict(g1), predict(g2))
  expect_equal(predict(g1), scores(g1, newdata=nuclearplants))
  expect_equal(predict(g2), scores(g2, newdata=nuclearplants))
  expect_equal(scores(g1, newdata=nuclearplants),
               scores(g2, newdata=nuclearplants))

  expect_equal(predict(g1), predict(g3))
  expect_equal(predict(g1), scores(g1, newdata=nuclearplants))
  expect_equal(predict(g3), scores(g3, newdata=nuclearplants))
  expect_equal(scores(g1, newdata=nuclearplants),
               scores(g3, newdata=nuclearplants))


})
