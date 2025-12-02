# (Internal) Predict for CBPS objects

The CBPS package fits ‘covariate balancing propensity score’ for use in
propensity score weighting. In the binary treatment case they can also
be used for matching. This method helps to implement that process in a
manner consistent with use of propensity scores elsewhere in optmatch;
see
[`scores`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
documentation.

## Usage

``` r
# S3 method for class 'CBPS'
predict(object, newdata = NULL, type = c("link", "response"), ...)
```

## Arguments

- object:

  A CBPS object

- newdata:

  Unused.

- type:

  Return inverse logits of fitted values (the default) or fitted values
  themselves

- ...:

  Unused.

## Value

Inverse logit of the fitted values.
