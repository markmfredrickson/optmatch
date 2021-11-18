test_that("strata function", {
  s <- strata(1:10, 1:10)
  expect_is(s, "character")
  expect_length(s, 10)
  expect_length(unique(s), 10)

  s <- strata(c(1, 1, 2, 2), c(1, 2, 2, 2))
  expect_length(s, 4)
  expect_length(unique(s), 3)

  s <- strata(c(1, NA, 1, NA), c(1, 1, NA, NA))
  expect_equal(is.na(s), c(FALSE, TRUE, TRUE, TRUE))

  expect_error(strata(1, 1:2), "same length")
})
