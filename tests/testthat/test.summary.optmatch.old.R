
context('summary optmatch')

test_that("summary.optmatch", {
  data(plantdist)
  expect_warning(s1 <- summary(f1 <- fullmatch(1 * (plantdist < 10)))) # a zero-1 matrix
  expect_true(all.equal(s1$thematch, f1))
  expect_true(is.null(s1$matching.failed))
  expect_true(all.equal(as.vector(s1$matched.set.structures), c(5,1,1)))
  expect_equal(s1$effective.sample.size, 8.1794871)
  #expect_equal(s1$total.distance, 0)
  #expect_equal(s1$total.tolerances, .0054166666)
  #expect_equal(sum(s1$matched.dist.quantiles), 0)


  # Mtching doesn't fail everywhere
  #expect_error(summary(pairmatch(plantdist + caliper(plantdist, 1)))) # Matching fails everywhere

  data(nuclearplants)
  psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
  psd <- match_on(psm, standardization.scale = sd) # backwards compatible to 0.7-2
  psfm <- fullmatch(psd + caliper(psd, 0.25), data = nuclearplants)
  summary(psfm) #!

  # Matching fails in a subgroup
  pspm <- pairmatch(caliper(match_on(psm, standarization.scale = sd,
                                     within = exactMatch(pr ~ pt, data = nuclearplants)),
                            width=2),
                    data=nuclearplants)

  expect_true(!is.null(summary(pspm)$matching.failed))
  psd[1,] <- psd[1,] + rep(100,22)

  # due to slight differences in the match on different platforms, just check that the
  # total.distances remain the same
  #expect_equal(summary(pairmatch(psd, controls=2, data = nuclearplants))$total.distance, 225.83338)

  # RItools is loaded directly, so this occasion can not happen
  # without PEBKAC.
  ## if ("RItools" %in% loadedNamespaces()) {
  ##   detach(package:RItools, unload=TRUE)
  ## }
  ## s2 <- summary(psfm, propensity.model=psm)
  ## expect_true(!is.null(s2$warnings))

  require('RItools')
  s3 <- summary(psfm, propensity.model='foo')
  expect_true(!is.null(s3$warnings))
  s4 <- summary(psfm, propensity.model=psm)
  expect_true(is.null(s4$warnings))
  s5 <- summary(psfm, psm)
  expect_true(is.null(s5$warnings))

  #expect_equal(s2$thematch, s3$thematch)
  #expect_equal(s2$thematch, s4$thematch)
  #expect_equal(s2$thematch, s5$thematch)


  psm2 <- glm(pr~ cut(date, c(67, 69.5, 72)) +
              t1 + t2 + cap + ne + ct + bw + cum.n + pt,
              family=binomial, data=nuclearplants)
  psd2 <- match_on(psm2, standardization.scale = sd)
  psd2summary <- summary(pairmatch(psd2, data = nuclearplants), propensity.model=psm2)

  # due to slight differences in the match on different platforms, just check that the
  # total.distances are the same and that the chi-squared value is 9.5 +- 0.5

  #expect_equal(psd2summary$total.distance, 7.5621504)
  chisquared.value <- psd2summary$balance$overall$chisquare
  expect_true(abs(9.5 - chisquared.value) < 0.5)
})
