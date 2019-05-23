################################################################################
## Tests to ensure that results are in fact optimal
################################################################################

context("optimality")

test_that("Simple example with and without calipers", {
    x <- data.frame(row.names = c("A", "B", "C", "D", "E"), z = c(1, 1, 0, 0, 0), y = c(0, 2, -3, 0, 6))
    m <- match_on(z ~ y, data = x, method = "euclidean")

    pm <- pairmatch(m, data = x)
    expect_true(pm["A"] == pm["D"])
    expect_true(pm["B"] == pm["E"])

    pmc4 <- pairmatch(caliper(m, 4, values = TRUE), data = x)
    expect_equivalent(pmc4, pm)

    fm <- fullmatch(m, data = x)
    fmn <- as.numeric(fm)
    expect_true(fm["A"] == fm["C"] && fm["A"] == fm["D"])
    expect_true(fm["B"] == fm["E"] && fm["B"] != fm["A"])
})
