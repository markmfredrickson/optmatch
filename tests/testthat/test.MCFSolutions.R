context("MCFSolutions & co: S4 classes to encode min-cost-flow solutions")

test_that("Instantiation & validity", {
    mbls1  <- new("MatchablesInfo", data.frame(name=c('a','b'), subproblem=rep("c",2),
                                              'kind'=c('treatment', 'control'), stringsAsFactors=F))
    expect_true(validObject(mbls1))

    mbls2  <- mbls1
    colnames(mbls2)[3]  <- "Kind"
    expect_error(validObject(mbls2), "eed a 'kind' column", fixed=TRUE)

    mbls3  <- mbls1
    mbls3@.Data[[3]]  <- rep("Exposed", 2)
    expect_error(validObject(mbls3), "values other than 'treatment' or 'control'", fixed=TRUE)
})
