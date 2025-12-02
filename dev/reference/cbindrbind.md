# Combine InfinitySparseMatrices or BlockedInfinitySparseMatrices by row or column

This matches the syntax and semantics of cbind and rbind for matrices.

## Usage

``` r
# S3 method for class 'InfinitySparseMatrix'
cbind(x, y, ...)

# S3 method for class 'InfinitySparseMatrix'
rbind(x, y, ...)

# S3 method for class 'BlockedInfinitySparseMatrix'
cbind(x, y, ...)

# S3 method for class 'BlockedInfinitySparseMatrix'
rbind(x, y, ...)
```

## Arguments

- x:

  An InfinitySparseMatrix or BlockedInfinitySparseMatrix, agreeing with
  `y` in the appropriate dimension.

- y:

  An InfinitySparseMatrix or BlockedInfinitySparseMatrix, agreeing with
  `x` in the appropriate dimension.

- ...:

  Other arguments ignored.

## Value

A combined InfinitySparseMatrix or BlockedInfinitySparseMatrix

## Author

Mark Fredrickson
