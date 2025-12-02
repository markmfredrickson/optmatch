# Compute the effective sample size of a match.

The effective sample size is the sum of the harmonic means of the number
units in treatment and control for each matched group. For k matched
pairs, the effective sample size is k. As matched groups become more
unbalanced, the effective sample size decreases.

## Usage

``` r
effectiveSampleSize(x, z = NULL)

# S3 method for class 'factor'
effectiveSampleSize(x, z = NULL)

# Default S3 method
effectiveSampleSize(x, z = NULL)

# S3 method for class 'table'
effectiveSampleSize(x, z = NULL)
```

## Arguments

- x:

  An `optmatch` object, the result of
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

- z:

  A treatment indicator, a vector the same length as `match`. This is
  only required if the `match` object does not contain the
  contrast.group' attribute.

## Value

The equivalent number of pairs in this match.

## See also

[`summary.optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md),
[`stratumStructure`](https://markmfredrickson.github.io/optmatch/dev/reference/stratumStructure.md)
