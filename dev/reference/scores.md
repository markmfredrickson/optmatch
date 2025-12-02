# Extract scores (propensity, prognostic,...) from a fitted model

This is a wrapper for `predict`, adapted for use in matching. Given a
fitted model but no explicit `newdata` to ‘predict’ from, it constructs
its own `newdata` in a manner that's generally better suited for
matching.

## Usage

``` r
scores(object, newdata = NULL, ...)
```

## Arguments

- object:

  fitted model object determining scores to be generated.

- newdata:

  (optional) data frame containing variables with which scores are
  produced.

- ...:

  additional arguments passed to `predict`.

## Value

See individual `predict` functions.

## Details

Like `predict`, its default predictions from a `glm` are on the scale of
the linear predictor, not the scale of the response; see Rosenbaum \\
Rubin (1985). (This default can be overridden by specifying
`type="response"`.) In contrast to `predict`, if `scores` isn't given an
explicit `newdata` argument then it attempts to reconstruct one from the
context in which it is called, rather than from its first argument. For
example, if it's called within the `formula` argument of a call to
`glm`, its `newdata` is the same data frame that `glm` evaluates that
formula in, as opposed to the model frame associated with `object`. See
Examples.

The handling of missing independent variables also differs from that of
`predict` in two ways. First, if the data used to generate `object` has
`NA` values, they're mean-imputed using
[`fill.NAs`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md).
Secondly, if `newdata` (either the explicit argument, or the implicit
data generated from `object`) has `NA` values, they're likewise
mean-imputed using
[`fill.NAs`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md).
Also, missingness flags are added to the formula of `object`, which is
then re-fit, using
[`fill.NAs`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md),
prior to calling `predict`.

If `newdata` is specified and contains no missing data, `scores` returns
the same value as `predict`.

## References

P.~R. Rosenbaum and D.~B. Rubin (1985), ‘Constructing a control group
using multivariate matched sampling methods that incorporate the
propensity score’, *The American Statistician*, **39** 33–38.

## See also

[`predict`](https://rdrr.io/r/stats/predict.html)

## Author

Josh Errickson

## Examples

``` r
data(nuclearplants)
pg <- lm(cost~., data=nuclearplants, subset=(pr==0))
# The following two lines produce identical results.
ps1 <- glm(pr~cap+date+t1+bw+predict(pg, newdata=nuclearplants),
           data=nuclearplants)
#> Warning: prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
ps2 <- glm(pr~cap+date+t1+bw+scores(pg), data=nuclearplants)
#> Warning: prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
```
