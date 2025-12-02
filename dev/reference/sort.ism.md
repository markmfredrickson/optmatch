# Sort the internal structure of an InfinitySparseMatrix.

Internally, an `InfinitySparseMatrix` (Blocked or non) comprises of
vectors of values, row positions, and column positions. The ordering of
these vectors is not enforced. This function sorts the internal
structure, leaving the external structure unchanged (e.g.
`as.matrix(ism)` and `as.matrix(sort(ism))` will look identical despite
sorting.)

## Usage

``` r
# S3 method for class 'InfinitySparseMatrix'
sort(x, decreasing = FALSE, ..., byCol = FALSE)

# S3 method for class 'BlockedInfinitySparseMatrix'
sort(x, decreasing = FALSE, ..., byCol = FALSE)
```

## Arguments

- x:

  An `InfinitySparseMatrix` or `BlockedInfinitySparseMatrix`.

- decreasing:

  Logical. Should the sort be increasing or decreasing? Default `FALSE`.

- ...:

  Additional arguments ignored.

- byCol:

  Logical. Defaults to `FALSE`, so the returned ISM is row-dominant.
  `TRUE` returns a column-dominant ISM.

## Value

An object of the same class as `x` which is sorted according to `byCol`.

## Details

By default, the `InfinitySparseMatrix` is row-dominant, meaning the row
positions are sorted first, then column positions are sorted within each
row. Use argument `byCol` to change this.
