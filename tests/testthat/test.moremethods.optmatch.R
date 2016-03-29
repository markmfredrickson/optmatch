
context("more optmatch methods")

test_that("", {
  data(plantdist)
  expect_warning(plantsfm <- fullmatch(plantdist))
  p1 <- plantsfm[1:10]
  #a1 <- attributes(plantsfm[1:10])
  expect_equal(names(p1), c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))
  expect_equal(levels(p1), c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7"))
  expect_true("optmatch" %in% class(p1))
  expect_equal(attr(p1,"contrast.group"), c(rep(TRUE, 7), rep(FALSE, 3)))
  #expect_equal(unname(unlist(attr(p1, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))

  p2 <- plantsfm[5:10]
  expect_equal(names(p2), c("E", "F", "G", "H", "I", "J"))
  expect_equal(levels(p2), c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7"))
  expect_true("optmatch" %in% class(p2))
  expect_equal(attr(p2,"contrast.group"), c(rep(TRUE, 3), rep(FALSE, 3)))
  #expect_equal(unname(unlist(attr(p2, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))

  expect_true(all.equal(plantsfm[1:26 < 11], plantsfm[1:10]))

  p3 <- plantsfm[1:26 <6]
  expect_equal(names(p3), c("A", "B", "C", "D", "E"))
  expect_equal(levels(p3), c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7"))
  expect_true("optmatch" %in% class(p3))
  expect_equal(attr(p3,"contrast.group"), rep(TRUE,5))
  #expect_equal(unname(unlist(attr(p3, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))

  p4 <- plantsfm[1:26 <6, drop=TRUE]
  expect_equal(p3, p4, check.attributes=FALSE)
  expect_equal(levels(p4), c("1.1", "1.2", "1.3", "1.4", "1.5"))
  expect_equal(attributes(p3)[c(1,3,4)], attributes(p4)[c(1,3,4)])

  expect_true(all.equal(plantsfm[1:26 < 11], plantsfm[names(plantsfm)[1:10] ]))

  expect_equal(attributes(p1), attributes(plantsfm[names(plantsfm)[1:10] ]))

  p5 <- plantsfm[names(plantsfm)[5:10],drop=TRUE]
  expect_equal(p2, p5, check.attributes=FALSE)
  expect_equal(levels(p5), c("1.1", "1.2", "1.4", "1.5", "1.6", "1.7"))
  expect_equal(attributes(p2)[c(1,3,4)], attributes(p5)[c(1,3,4)])


  plantsfm[5] <- "1.4"
  p6 <- plantsfm[1:5]
  expect_equal(attributes(p3), attributes(p6))
  expect_true(!all(p3 == p6))

  expect_warning(plantsfm <- fullmatch(plantdist))
  p7 <- plantsfm
  p8 <- plantsfm[26:1]
  expect_equal(names(p7)[26:1], names(p8))
  expect_equal(levels(p7), levels(p8))


  ### arises in lme4:::lmerFactorList , which is called in lme4::lmer
  ### at following line:
  ###
  ###   fl <- lapply(bars, function(x) eval(substitute(as.factor(fac)[,
  ###        drop = TRUE], list(fac = x[[3]])), mf))
  ###
  ### (caused [.optmatch to die in optmatch version 0.4-0 on R-2.6.0 +)
  p9 <- as.factor(plantsfm)[, drop = TRUE]
  p10 <- plantsfm[, drop = TRUE]
  p11 <- plantsfm[, drop = FALSE]
  expect_true(all(p9==p10))
  expect_true(all(p9==p11))
})
