# optmatch_restrictions

Returns the restrictions which were used to generate the match.

## Usage

``` r
optmatch_restrictions(obj)
```

## Arguments

- obj:

  An optmatch object

## Value

A list of `min.controls`, `max.controls` and either `omit.fraction` or
`mean.controls`.

## Details

If `mean.controls` was explicitly specified in the creation of the
optmatch object, it is returned; otherwise `omit.fraction` is given.

Note that if the matching algorithm attempted to recover from initial
infeasible restrictions, the output from this function may not be the
same as the original function call.
