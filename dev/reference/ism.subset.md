# Subsetting for InfinitySparseMatrices

This matches the syntax and semantics of subset for matrices.

## Usage

``` r
# S3 method for class 'InfinitySparseMatrix'
subset(x, subset, select, ...)

# S4 method for class 'InfinitySparseMatrix'
x[i, j = NULL, ..., drop = TRUE]

# S4 method for class 'InfinitySparseMatrix'
x[i, j] <- value
```

## Arguments

- x:

  InfinitySparseMatrix to be subset or bound.

- subset:

  Logical expression indicating rows to keep.

- select:

  Logical expression indicating columns to keep.

- ...:

  Other arguments are ignored.

- i:

  Row indices.

- j:

  Col indices.

- drop:

  Ignored.

- value:

  replacement values

## Value

An InfinitySparseMatrix with only the selected elements.

## Author

Mark Fredrickson
