# (Deprecated, in favor of [`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)) Create matching distances

Deprecated in favor of
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)

## Usage

``` r
mdist(x, structure.fmla = NULL, ...)

# S3 method for class 'optmatch.dlist'
mdist(x, structure.fmla = NULL, ...)

# S3 method for class '`function`'
mdist(x, structure.fmla = NULL, data = NULL, ...)

# S3 method for class 'formula'
mdist(x, structure.fmla = NULL, data = NULL, subset = NULL, ...)

# S3 method for class 'glm'
mdist(x, structure.fmla = NULL, standardization.scale = mad, ...)

# S3 method for class 'bigglm'
mdist(x, structure.fmla = NULL, data = NULL, standardization.scale = mad, ...)

# S3 method for class 'numeric'
mdist(x, structure.fmla = NULL, trtgrp = NULL, ...)
```

## Arguments

- x:

  The object to use as the basis for forming the mdist. Methods exist
  for formulas, functions, and generalized linear models.

- structure.fmla:

  A formula denoting the treatment variable on the left hand side and an
  optional grouping expression on the right hand side. For example,
  `z ~ 1` indicates no grouping. `z ~ s` subsets the data only computing
  distances within the subsets formed by `s`. See method notes, below,
  for additional formula options.

- ...:

  Additional method arguments. Most methods require a 'data' argument.

- data:

  Data where the variables references in `x` live.

- subset:

  If non-`NULL`, the subset of `data` to be used.

- standardization.scale:

  A function to scale the distances; by default uses `mad`.

- trtgrp:

  Dummy variable for treatment group membership.

## Value

Object of class `optmatch.dlist`, which is suitable to be given as
`distance` argument to
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
or
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

## Details

The `mdist` method provides three ways to construct a matching distance
(i.e., a distance matrix or suitably organized list of such matrices):
guided by a function, by a fitted model, or by a formula. The class of
the first argument given to `mdist` determines which of these methods is
invoked.

The `mdist.function` method takes a function of two arguments. When
called, this function will receive the treatment observations as the
first argument and the control observations as the second argument. As
an example, the following computes the raw differences between values of
`t1` for treatment units (here, nuclear plants with `pr==1`) and
controls (here, plants with `pr==0`), returning the result as a distance
matrix:

`` sdiffs <- function(treatments, controls) { abs(outer(treatments$t1, controls$t1, `-`)) }  ``

The `mdist.function` method does similar things as the earlier optmatch
function `makedist`, although the interface is a bit different.

The `mdist.formula` method computes the squared Mahalanobis distance
between observations, with the right-hand side of the formula
determining which variables contribute to the Mahalanobis distance. If
matching is to be done within strata, the stratification can be
communicated using either the `structure.fmla` argument (e.g. `~ grp`)
or as part of the main formula (e.g. `z ~ x1 + x2 | grp`).

An `mdist.glm` method takes an argument of class `glm` as first
argument. It assumes that this object is a fitted propensity model,
extracting distances on the linear propensity score (logits of the
estimated conditional probabilities) and, by default, rescaling the
distances by the reciprocal of the pooled s.d. of treatment- and
control-group propensity scores. (The scaling uses `mad`, for resistance
to outliers, by default; this can be changed to the actual s.d., or
rescaling can be skipped entirely, by setting argument
`standardization.scale` to `sd` or `NULL`, respectively.) A
`mdist.bigglm` method works analogously with `bigglm` objects, created
by the `bigglm` function from package ‘biglm’, which can handle bigger
data sets than the ordinary glm function can. In contrast with
`mdist.glm` it requires additional `data` and `structure.fmla`
arguments. (If you have enough data to have to use `bigglm`, then you'll
probably have to subgroup before matching to avoid memory problems. So
you'll have to use the `structure.fmla` argument anyway.)

## References

P.~R. Rosenbaum and D.~B. Rubin (1985), ‘Constructing a control group
using multivariate matched sampling methods that incorporate the
propensity score’, *The American Statistician*, **39** 33–38.

## See also

[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md),
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)

## Author

Mark M. Fredrickson
