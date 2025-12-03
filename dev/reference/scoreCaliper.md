# (Internal) Helper function to create an InfinitySparseMatrix from a set of scores, a treatment indicator, and a caliper width.

(Internal) Helper function to create an InfinitySparseMatrix from a set
of scores, a treatment indicator, and a caliper width.

## Usage

``` r
scoreCaliper(x, z, caliper, within = NULL)
```

## Arguments

- x:

  The scores, a vector indicating the 1-D location of each unit.

- z:

  The treatment assignment vector (same length as `x`)

- caliper:

  The width of the caliper with respect to the scores `x`.

- within:

  A valid distance specification, such as the result of
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
  or
  [`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md).
  Finite entries indicate which distances to create. Including this
  argument can significantly speed up computation for sparse matching
  problems. Specify this filter either via `within` or via `strata`
  elements of a formula; mixing these methods will fail.

## Value

An `InfinitySparseMatrix` object, suitable to be passed to
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
as an `within` argument.
