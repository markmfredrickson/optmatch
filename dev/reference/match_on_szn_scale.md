# pooled dispersion for a numeric variable

Pooled Dispersion for a Numeric Variable

## Usage

``` r
match_on_szn_scale(x, trtgrp, standardizer = mad)
```

## Arguments

- x:

  numeric variable

- trtgrp:

  logical or numeric. If numeric, coerced to logical via `!`

- standardizer:

  function or numeric of length 1

## Value

numeric of length 1

## Details

Dispersion as pooled across a treatment and a control group. By default,
the measure of dispersion calculated within each group is not the
ordinary standard deviation but rather the robust alternative provided
by [`stats::mad`](https://rdrr.io/r/stats/mad.html).
