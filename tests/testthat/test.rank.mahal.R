################################################################################
# compute_rank.mahalanobis: calculate the mahalanobis distance between ranks
# of vectors
################################################################################

context("compute_rank.mahalanobis function")

if (requireNamespace("MASS")) {
# Rosenbaum's original rank mahalanobis distance code
# from Design of observational studies
# http://www-stat.wharton.upenn.edu/~rosenbap/Rdospublic.RData
# Downloaded 2013.06.25
# Subsequently touched up a tad, by @benthestatistician
      compute_smahal <- function (z, X) {
        X <- as.matrix(X)
        n <- dim(X)[1]
        rownames(X) <- 1:n
        k <- dim(X)[2]
        m <- sum(z)
        for (j in 1:k) X[, j] <- rank(X[, j])
        cv <- cov(X)
        vuntied <- var(1:n)
        rat <- sqrt(vuntied/diag(cv))
        cv <- diag(rat) %*% cv %*% diag(rat)
        out <- matrix(NA, m, n - m)
        Xc <- X[z == 0, , drop=FALSE]
        Xt <- X[z == 1, , drop=FALSE]
        rownames(out) <- rownames(X)[z == 1]
        colnames(out) <- rownames(X)[z == 0]
        icov <- MASS::ginv(cv)
        for (i in 1:m) out[i, ] <- mahalanobis(Xc, Xt[i, ], icov,
                                               inverted = T)
        out
    }

}

test_that('compute_rank.mahal returns results similar to the Rosenbaum code', {

  if (requireNamespace("MASS")) {
    nr <- 10L
    nc <- 5L

    z <- integer(nr)
    z[sample(1:nr, nr / 2L)] <- 1L
    numdists <- sum(z)*sum(!z)
      
    X <- matrix(runif(nr * nc), nr)
    row.names(X) <- 1L:nr
      
    reference_rankmahal <- compute_smahal(z, X)
    expect_that(optmatch:::compute_rank_mahalanobis(NULL, X, as.logical(z)),
                is_equivalent_to(reference_rankmahal))

      
    expect_equivalent(as.matrix(match_on(z ~ X, method = "rank")), reference_rankmahal)

    df <- data.frame(z = z, X)
    expect_equivalent(as.matrix(match_on(z ~ ., data = df, method = "rank")), reference_rankmahal)
  }
})

test_that("Fix for #128 (`compute_rank_mahalanobis` ignores index argument) holds", {

      if (requireNamespace("MASS")) {
    nr <- 10L
    nc <- 5L

    z <- integer(nr)
    z[sample(1:nr, nr / 2L)] <- 1L
    numdists <- sum(z)*sum(!z)
      
    X <- matrix(runif(nr * nc), nr)
    row.names(X) <- 1L:nr
      
    reference_rankmahal <- compute_smahal(z, X)

        indices <- expand.grid(rownames(reference_rankmahal), colnames(reference_rankmahal))
    indices <- as.matrix(indices)
    expect_equivalent(optmatch:::compute_rank_mahalanobis(indices, X, as.logical(z)),
                reference_rankmahal[1L:numdists])

    expect_equivalent(optmatch:::compute_rank_mahalanobis(indices[-1L, ], X, as.logical(z)),
                      reference_rankmahal[-1L])
    omitdists <- sample(1L:numdists, 2)
    expect_equivalent(optmatch:::compute_rank_mahalanobis(indices[-omitdists, ],
                                                          X, as.logical(z)),
                      reference_rankmahal[-omitdists])

      }
})
