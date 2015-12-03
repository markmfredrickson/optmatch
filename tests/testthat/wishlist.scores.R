#############################################################################
# wishlist for scores
#############################################################################

library(testthat)

context("scores wishlist")

test_that("scores works with attach/detach", {
  pg <- lm(cost ~ ., data=nuclearplants, subset=(pr==0))
### function as predict correctly
  pred <- predict(pg, newdata=nuclearplants)
  attach(nuclearplants)
  scores <- scores(pg)
  detach(nuclearplants)
  expect_equal(pred, scores, check.attributes=FALSE)

})

test_that("scores works with attach/detach in a model", {
  pg <- lm(cost~., data=nuclearplants, subset=(pr==0))

  ps1 <- glm(pr ~ cap + date + t1 + bw + predict(pg, newdata=nuclearplants), data=nuclearplants)
  ps2 <- with(nuclearplants, glm(pr ~ cap + date + t1 + bw + scores(pg)))
  attach(nuclearplants)
  expect_error(ps3 <- glm(pr ~ cap + date + t1 + bw + predict(pg)))
  ps4<- glm(pr ~ cap + date + t1 + bw + scores(pg))
  detach(nuclearplants)

  expect_equal(fitted(ps1), fitted(ps2), check.attributes=FALSE)
  expect_equal(fitted(ps1), fitted(ps4), check.attributes=FALSE)
})
