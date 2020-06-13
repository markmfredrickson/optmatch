context("edgelist functions")

test_that("Basic", {
    cns <- c("C1", "C2", "C3", "C4")
    rns <- c("TZ", "TY", "TX", "TW", "TU")
    ism <- makeInfinitySparseMatrix(1:7,
                                    cols = c(1, 2, 1, 2, 3, 1, 3),
                                    rows = c(2, 2, 3, 3, 3, 5, 5),
                                    colnames = cns,
                                    rownames = rns)
    ## ISM method
    ism.el <- edgelist(ism, c(rns, cns))

    expect_is(ism.el, "EdgeList")
    expect_equal(subset(ism.el, i == "TU" & j == "C3")$dist, 7)

    ## matrix method
    m <- as.matrix(ism)
    m.el <- edgelist(m, c(rns, cns))

    expect_is(m.el, "EdgeList")
    expect_equal(subset(m.el, i == "TU" & j == "C3")$dist, 7)
    ## EdgeList method
    expect_identical(edgelist(m.el, c(rns, cns)), m.el)

    ## data frame method (i.e. an EdgeList that had its S4 class dropped)
    expect_false(is(subset(m.el, TRUE), "EdgeList"))
    expect_is(subset(m.el, TRUE), "data.frame")
    expect_identical(edgelist(subset(m.el, TRUE)), m.el)

    ## error on other classes
    expect_error(edgelist(numeric(1), c(rns, cns)), "Not implemented")
})
test_that("default levels",{
    cns <- c("C1", "C2", "C3", "C4")
    rns <- c("TZ", "TY", "TX", "TW", "TU")
    ism <- makeInfinitySparseMatrix(1:7,
                                    cols = c(1, 2, 1, 2, 3, 1, 3),
                                    rows = c(2, 2, 3, 3, 3, 5, 5),
                                    colnames = cns,
                                    rownames = rns)
    ism.el <- edgelist(ism, c(rns, cns))
    expect_identical(edgelist(ism), ism.el)


    m <- as.matrix(ism)
    m.el <- edgelist(m, c(rns, cns))
    expect_identical(edgelist(m), m.el)
    ## (too beat to test this for the other methods)
})
test_that("remove edges to or from nodes not appearing in y",{
    cns <- c("C1", "C2", "C3", "C4")
    rns <- c("TZ", "TY", "TX", "TW", "TU")
    ism <- makeInfinitySparseMatrix(1:7,
                                    cols = c(1, 2, 1, 2, 3, 1, 3),
                                    rows = c(2, 2, 3, 3, 3, 5, 5),
                                    colnames = cns,
                                    rownames = rns)
    ## ISM method
    ism.el <- edgelist(ism, c(rns[-2], cns))
    expect_equal(nrow(ism.el), length(ism)-2)
    ## verify no scrambling of values:
    expect_equal(subset(ism.el, i == "TU" & j == "C3")$dist, 7)

    ## matrix method
    m <- as.matrix(ism)
    m.el <- edgelist(m, c(rns[-2], cns))
    expect_equal(nrow(m.el), length(ism)-2)
    ## verify no scrambling of values:
    expect_equal(subset(m.el, i == "TU" & j == "C3")$dist, 7)

    expect_is(m.el, "EdgeList")
    expect_equal(subset(m.el, i == "TU" & j == "C3")$dist, 7)
    ## EdgeList method
    expect_identical(edgelist(m.el, c(rns[-2], cns)), m.el)

    ## data frame method (i.e. an EdgeList that had its S4 class dropped)
    expect_false(is(subset(m.el, TRUE), "EdgeList"))
    expect_is(subset(m.el, TRUE), "data.frame")
    expect_equivalent(edgelist(subset(m.el, TRUE), c(rns[-2], cns)), m.el)
})
