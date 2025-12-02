# pooled dispersion for a numeric variable

Dispersion as pooled across a treatment and a control group. By default,
the measure of dispersion calculated within each group is not the
ordinary standard deviation as in
[`stats::sd`](https://rdrr.io/r/stats/sd.html) but rather the robust
alternative encoded in [`stats::mad`](https://rdrr.io/r/stats/mad.html).
The dispersion measurements are combined by squaring, averaging with
weights proportional to one minus the sizes of the groups and then
taking square roots. Used in
[`match_on.glm`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md).

## Usage

``` r
standardization_scale(x, trtgrp, standardizer = NULL, svydesign_ = NULL)
```

## Arguments

- x:

  numeric variable

- trtgrp:

  logical or numeric. If numeric, coerced to logical via `!`

- standardizer:

  function, `NULL` or numeric of length 1

- svydesign\_:

  ordinarily `NULL`, but may also be a `survey.design2`; see Details.

## Value

numeric of length 1

## Details

A non-NULL `svydesign_` parameter indicates that the dispersion
calculations are to be made respecting the weighting scheme implicit in
that `survey.design2` object. If `standardizer` is `NULL`, one gets a
calculation in the style of
[`stats::mad`](https://rdrr.io/r/stats/mad.html) but with weights,
performed by `optmatch:::svy_sd`; for a pooling of weighted standard
deviations, one would pass a non-`NULL` `svydesign_` parameter along
with `standardizer=optmatch:::svy_sd`. (More generally, the provided
`standardizer` function should accept as a sole argument a
`survey.design2` object, with `nrows(svydesign_$variables)` equal to the
lengths of `x` and `trtgrp`. This object is expected to carry a numeric
variable ‘`x`’, and the `standardizer` function is to return the
dispersion of this variable.)
