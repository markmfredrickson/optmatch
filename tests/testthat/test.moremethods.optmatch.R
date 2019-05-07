
context("more optmatch methods")
#CHANGING BEHAVIOR OF DROP ARGUMENT, SEE OptmatchS4.R @ line 153. no drop is the same as drop = TRUE under new expectations
test_that("", {
  data(plantdist)
  expect_warning(plantsfm <- fullmatch(plantdist))
  p1 <- plantsfm[1:10]
  #a1 <- attributes(plantsfm[1:10])
  expect_equal(names(p1), c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))
  expect_equal(levels(p1), c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7"))
  expect_true("Optmatch" %in% class(p1))
  expect_equal(p1@node.data[match(names(p1), p1@node.data$name), c("contrast.group")], c(rep(TRUE, 7), rep(FALSE, 3)))
  #expect_equal(unname(unlist(attr(p1, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))

  p2 <- plantsfm[5:10]
  expect_equal(names(p2), c("E", "F", "G", "H", "I", "J"))
  expect_equal(levels(p2), c("1.1", "1.2", "1.4", "1.5", "1.6", "1.7"))
  expect_true("Optmatch" %in% class(p2))
  expect_equal(p2@node.data[match(names(p2), p2@node.data$name), c("contrast.group")], c(rep(TRUE, 3), rep(FALSE, 3)))
  #expect_equal(unname(unlist(attr(p2, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))

  expect_true(all.equal(plantsfm[1:26 < 11]@.Data, plantsfm[1:10]@.Data))

  p3 <- plantsfm[1:26 <6]
  expect_equal(names(p3), c("A", "B", "C", "D", "E"))
  expect_equal(levels(p3), c("1.1", "1.2", "1.3", "1.4", "1.5"))
  expect_true("Optmatch" %in% class(p3))
  expect_equal(p3@node.data[match(names(p3), p3@node.data$name), c("contrast.group")], rep(TRUE,5))
  #expect_equal(unname(unlist(attr(p3, "matched.distances"))), c(0, 0, 4, 6, 9, 7, 8, 2, 6, 0, 4, 0, 2, 8, 4, 5, 12, 4, 8))


  p4 <- plantsfm[1:26 <6, drop=TRUE]
  expect_equal(p3@.Data, p4@.Data, check.attributes=FALSE)
  expect_equal(levels(p4), c("1.1", "1.2", "1.3", "1.4", "1.5"))
  expect_equal(attributes(p3)[c(3,5)], attributes(p4)[c(3,5)])
  expect_true(class(p3) == "Optmatch")
  expect_true(class(p4) == "Optmatch")

  all.equal(as.optmatch(plantsfm[1:26 < 11]), as.optmatch(plantsfm[names(plantsfm)[1:10] ]))

  expect_equal(attributes(p1), attributes(plantsfm[names(plantsfm)[1:10] ]))

  p5 <- plantsfm[names(plantsfm)[5:10],drop=TRUE]
  #NOTE: I AM COMMENTING OUT THE FOLLOWING TEST BECAUSE I DON'T BELIEVE IT SHOULD NECESSARILY BE EXPECTED TO BE TRUE UNDER THE NEW DATA STRUCTURE ARRANGEMENT
  #expect_equal(p2@, p5@, check.attributes=FALSE)
  expect_equal(levels(p5), c("1.1", "1.2", "1.4", "1.5", "1.6", "1.7"))

  expect_equal(names(p2), names(p5))
  expect_equal(p2@node.data[match(names(p2), p2@node.data$name), c('contrast.group')], p5@node.data[match(names(p5), p5@node.data$name), c('contrast.group')])

  # Not sure there's a strong equivalent in the new structure for this kind of operation...
  t <- as.numeric(p3@.Data)
  t[5] <- t[5] + 2
  p6 <- plantsfm[1:5]
  expect_equal(attributes(p3), attributes(p6))
  expect_true(!all(t == p6@.Data))

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
