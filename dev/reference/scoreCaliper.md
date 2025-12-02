# (Internal) Helper function to create an InfinitySparseMatrix from a set of scores, a treatment indicator, and a caliper width.

(Internal) Helper function to create an InfinitySparseMatrix from a set
of scores, a treatment indicator, and a caliper width.

## Usage

``` r
scoreCaliper(x, z, caliper)
```

## Arguments

- x:

  The scores, a vector indicating the 1-D location of each unit.

- z:

  The treatment assignment vector (same length as `x`)

- caliper:

  The width of the caliper with respect to the scores `x`.

## Value

An `InfinitySparseMatrix` object, suitable to be passed to
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
as an `within` argument.
