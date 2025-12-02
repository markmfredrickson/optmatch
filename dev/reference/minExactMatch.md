# Find the minimal exact match factors that will be feasible for a given maximum problem size.

The
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
function creates a smaller matching problem by stratifying observations
into smaller groups. For a problem that is larger than maximum allowed
size, `minExactMatch` provides a way to find the smallest exact matching
problem that will allow for matching.

## Usage

``` r
minExactMatch(x, scores = NULL, width = NULL, maxarcs = 1e+07, ...)
```

## Arguments

- x:

  The object for dispatching.

- scores:

  Optional vector of scores that will be checked against a caliper
  width.

- width:

  Optional width of a caliper to place on the scores.

- maxarcs:

  The maximum problem size to attempt to fit.

- ...:

  Additional arguments for methods.

## Value

A factor grouping units, suitable for
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md).

## Details

`x` is a formula of the form `Z ~ X1 + X2`, where `Z` is indicates
treatment or control status, and `X1` and `X2` are variables can be
converted to factors. Any additional arguments are passed to
[`model.frame`](https://rdrr.io/r/stats/model.frame.html) (e.g., a
`data` argument containing `Z`, `X1`, and `X2`).

The the arguments `scores` and `width` must be passed together. The
function will apply the caliper implied by the scores and the width
while also adding in blocking factors.
