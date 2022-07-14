# test whether two matches are the same. Uses all.equal on exceedances
# to ignore errors below some tolerance. After checking those, strips
# attributes that may differ but not break `identical` status.
# prices are often not the name, so it is often useful to ignore those
# other arguments passed to expect_equal
match_equal <- function(match1, match2, ignore.solver = TRUE, ignore.prices = TRUE, ...) {
  expect_true(all.equal(attr(match1, "exceedances"),
                        attr(match2, "exceedances")))

  attr(match1, "hashed.distance") <- NULL
  attr(match2, "hashed.distance") <- NULL
  attr(match1, "exceedances") <- NULL
  attr(match2, "exceedances") <- NULL
  attr(match1, "call") <- NULL
  attr(match2, "call") <- NULL
  if (!ignore.solver) {
    attr(match1, "solver") <- NULL
    attr(match2, "solver") <- NULL
  }

  if (ignore.prices) {
    attr(match1, "MCFSolutions")@nodes$price <- NULL
    attr(match2, "MCFSolutions")@nodes$price <- NULL
  }

  expect_equal(match1, match2, ...)
}

#' Similar to match_equal, but doesn't care about differences
#' among labels of matched sets.
match_equivalent <- function(match1, match2) {
  expect_true(all.equal(attr(match1, "exceedances"),
                        attr(match2, "exceedances")))

  attr(match1, "hashed.distance") <- NULL
  attr(match2, "hashed.distance") <- NULL
  attr(match1, "exceedances") <- NULL
  attr(match2, "exceedances") <- NULL
  attr(match1, "call") <- NULL
  attr(match2, "call") <- NULL

  m1labs <- as.character(match1[!is.na(match1) & !duplicated(match1)])
  levels(match1)[match(m1labs, levels(match1))] <- m1labs
  match1 <- factor(match1)

  m2labs <- as.character(match2[!is.na(match2) & !duplicated(match2)])
  levels(match2)[match(m2labs, levels(match2))] <- m2labs
  match2 <- factor(match2)

  expect_true(identical(match1, match2))
}
