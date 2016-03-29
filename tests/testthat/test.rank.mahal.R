################################################################################
# compute_rank.mahalanobis: calculate the mahalanobis distance between ranks
# of vectors
################################################################################

context("compute_rank.mahalanobis function")

test_that('compute_rank.mahal returns results similar to the Rosenbaum code', {

  if (require(MASS)) {
    # Rosenbaum's original rank mahalanobis distance code
    # from Design of observational studies
    # http://www-stat.wharton.upenn.edu/~rosenbap/Rdospublic.RData
    # 2013.06.25
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
        Xc <- X[z == 0, ]
        Xt <- X[z == 1, ]
        rownames(out) <- rownames(X)[z == 1]
        colnames(out) <- rownames(X)[z == 0]
        icov <- ginv(cv)
        for (i in 1:m) out[i, ] <- mahalanobis(Xc, Xt[i, ], icov,
                                               inverted = T)
        out
    }

    nr <- 10L
    nc <- 5L

    z <- integer(nr)
    z[sample(1:nr, nr / 2L)] <- 1L

    X <- matrix(runif(nr * nc), nr)

    expect_that(compute_rank_mahalanobis(NULL, X, as.logical(z)),
                is_equivalent_to(compute_smahal(z, X)))

    expect_equivalent(as.matrix(match_on(z ~ X, method = "rank")), compute_smahal(z, X))

    df <- data.frame(z = z, X)
    expect_equivalent(as.matrix(match_on(z ~ ., data = df, method = "rank")), compute_smahal(z, X))
  }
})
