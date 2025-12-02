# Optmatch Class

The `optmatch` class describes the results of an optimal full matching
(using either
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
or
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)).
For the most part, these objects can be treated as `factors`.

The summary function quantifies `optmatch` objects on the effective
sample size, the distribution of distances between matched units, and
how well the match reduces average differences.

## Usage

``` r
# S3 method for class 'optmatch'
summary(
  object,
  propensity.model = NULL,
  ...,
  min.controls = 0.2,
  max.controls = 5,
  quantiles = c(0, 0.5, 0.95, 1)
)
```

## Arguments

- object:

  The `optmatch` object to summarize.

- propensity.model:

  An optional propensity model (the result of a call to `glm`) to use
  when summarizing the match. If the RItools package is installed, an
  additional chi-squared test will be performed on the average
  differences between treated and control units on each variable used in
  the model. See the `xBalance` function in the RItools package for more
  details.

- ...:

  Additional arguments to pass to `xBalance` when also passing a
  propensity model.

- min.controls:

  To minimize the the display of a groups with many treated and few
  controls, all groups with more than 5 treated units will be summarized
  as “5+”. This is the reciprocal of the default value (1/5 = 0.2).
  Lower this value to see more groups.

- max.controls:

  Like `min.controls` sets maximum group sized displayed with respect to
  the number of controls. Raise this value to see more groups.

- quantiles:

  A points in the ECDF at which the distances between units will be
  displayed.

## Value

`optmatch.summary`

## Details

`optmatch` objects descend from `factor`. Elements of this vector
correspond to members of the treatment and control groups in reference
to which the matching problem was posed, and are named accordingly; the
names are taken from the row and column names of `distance`. Each
element of the vector is either `NA`, indicating unavailability of any
suitable matches for that element, or the concatenation of: (i) a
character abbreviation of the name of the subclass (as encoded using
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md))
(ii) the string `.`; and (iii) a non-negative integer. In this last
place, positive whole numbers indicate placement of the unit into a
matched set and `NA` indicates that all or part of the matching problem
given to `fullmatch` was found to be infeasible. The functions
[`matched`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md),
[`unmatched`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md),
and
[`matchfailed`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md)
distinguish these scenarios.

Secondarily, `fullmatch` returns various data about the matching process
and its result, stored as attributes of the named vector which is its
primary output. In particular, the `exceedances` attribute gives upper
bounds, not necessarily sharp, for the amount by which the sum of
distances between matched units in the result of `fullmatch` exceeds the
least possible sum of distances between matched units in a feasible
solution to the matching problem given to `fullmatch`. (Such a bound is
also printed by `print.optmatch` and `summary.optmatch`.)

## See also

[`print.optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/print.optmatch.md)
