################################################################################
# Tests for makedist function
################################################################################

library(testthat)
context("compute_mahalanobis tests")

mahal.cov <- function(x, z) {
    mt <- cov(x[z, ,drop=FALSE]) * (sum(z) - 1) / (length(z) - 2)
    mc <- cov(x[!z, ,drop=FALSE]) * (sum(!z) - 1) / (length(!z) - 2)
    return(mt + mc)
}

trad.mahal <- function(index, x, z) {
    x.cov <- mahal.cov(x, z)
    ans <- apply(index, 1, function(pair)
        return(mahalanobis(x[pair[1], ], x[pair[2], ], x.cov)))
    return(sqrt(ans))
}

make.data <- function(k) {
    # generate some interesting but random vectors
    set.seed(20130811)

    x <- MASS::mvrnorm(k, mu = c(-1, 0, 1),
                  Sigma = matrix(c(1, 0.25, 0.25,
                      0.25, 1, 0.25,
                      0.25, 1, 0.25), nrow = 3))
    rownames(x) <- paste('v', 1:nrow(x), sep='')
    
    # top 10% assigned to treatment
    tmp <- rowSums(x)
    z <- vector("integer", k)
    z[order(tmp, decreasing = TRUE)[1:(0.1 * k)]] <- 1L
    z <- as.logical(z)
    
    # get treatment X control pairs into an index
    tns <- rownames(x)[z]
    nt <- length(tns)
    
    cns <- rownames(x)[!z]
    nc <- length(cns)
    
    t.coord <- rep(tns, nc)
    c.coord <- rep(cns, each = nt)

    index <- cbind(t.coord, c.coord)

    # return index, data, and z
    return(list(index = index, data = x, z = z))
}

test_that("Checking mahalanobis distance", {
    # random data set of size btw 100 and 500
    args <- make.data(sample(100:500, 1))
    
    expect_equal(
        compute_mahalanobis(args$index, args$data, args$z),
        trad.mahal(args$index, args$data, args$z))
})
