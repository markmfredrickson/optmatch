context("edgelist functions")

test_that("Basic", {
    library(dplyr)
    cns <- c("C1", "C2", "C3", "C4")
    rns <- c("TZ", "TY", "TX", "TW", "TU")
    ism <- makeInfinitySparseMatrix(1:7,
                                    cols = c(1, 2, 1, 2, 3, 1, 3),
                                    rows = c(2, 2, 3, 3, 3, 5, 5),
                                    colnames = cns,
                                    rownames = rns)
    ism.el <- edgelist(ism, c(rns, cns))

    expect_equal(filter(ism.el, i == "TU" & j == "C3")$dist, 7)


    m <- as.matrix(ism)
    m.el <- edgelist(m, c(rns, cns))

    expect_equal(filter(m.el, i == "TX" & j == "C3")$dist, 5)
})

test_that("remove edges to or from nodes not appearing in y",{
    library(dplyr)
    cns <- c("C1", "C2", "C3", "C4")
    rns <- c("TZ", "TY", "TX", "TW", "TU")
    ism <- makeInfinitySparseMatrix(1:7,
                                    cols = c(1, 2, 1, 2, 3, 1, 3),
                                    rows = c(2, 2, 3, 3, 3, 5, 5),
                                    colnames = cns,
                                    rownames = rns)
    ism.el <- edgelist(ism, c(rns[-2], cns))
    expect_equal(nrow(ism.el), length(ism)-2)
    expect_equal(filter(ism.el, i == "TU" & j == "C3")$dist, 7)
})
