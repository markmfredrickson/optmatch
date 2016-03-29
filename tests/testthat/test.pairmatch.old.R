
context("Pairmatch function old")

test_that("pairmatch", {
  data(plantdist)

  stripCall <- function(obj) {
    attr(obj, "call") <- NULL
    obj
  }


  expect_warning(p1 <- pairmatch(plantdist))
  expect_warning(f1 <- fullmatch(plantdist, max.controls = 1, min.controls = 1, omit.fraction = 1 - 7/19))
  expect_equal(stripCall(f1), stripCall(p1))

  expect_warning(p2 <- pairmatch(plantdist, controls=2))
  expect_warning(f2 <- fullmatch(plantdist, max.controls = 2, min.controls = 2, omit.fraction = 1 - 14/19))
  expect_equal(stripCall(f2), stripCall(p2))

  # plantdist has several 0's, so it won't fail everywhere.
#  expect_error(pairmatch(plantdist + caliper(plantdist, 1))) # Matching fails everywhere

  expect_warning(p3 <- pairmatch(plantdist + caliper(plantdist, 5, compare = `<`),
                                 remove.unmatchables=TRUE)) # Matching works after removing plant 'F'
  expect_warning(f3 <- fullmatch(plantdist + caliper(plantdist, 5, compare = `<`),
                                 max.controls = 1, min.controls = 0, omit.fraction = 1 - 6/19))

  expect_true(all.equal(f3,p3, check.attributes=FALSE))

  data(nuclearplants)
  # in both of match_on calls below use sd to maintain backwards compatibility with
  # pscore.dist, which used sd by default. match_on has used mad as the std. scale
  # since it was added to the package, so the use of match_on should be consistent
  # for users going forward.
  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  psd <- match_on(psm, standardization.scale = sd)
  pm <- pairmatch(psd, controls=2, data = nuclearplants)

  # the pm match immediately above was giving slightly different answers in some environment
  # the problem allowed multiple optimal solutions, and different choices were picked in different environments
  # the sum of matched distances should be the same across all environments

  # expect_true(all.equal(summary(pm)$total.distance, 25.83338, tolerance=1e-5))

  ## We no longer throw an error when a subclass fails
  # again an error would be thrown (which R CMD CHECK does not like)
  #expect_error(pairmatch(caliper(match_on(psm, standardization.scale = sd,
  #   within = exactMatch(pr ~ pt, data =
  #   nuclearplants)), width=2))) # Fails in subclass '1'
})
