# (Internal) A helper function to turn formulas into treatment and blocking variables

Given a function and any of the arguments normally passed to
model.frame, this function will return a data.frame with two columns: a
treatment indicator and a blocking factor.

## Usage

``` r
fmla2treatedblocking(x, ...)
```

## Arguments

- x:

  A formula

- ...:

  Arguments to be passed to model.frame (e.g. `data`)

## Value

data.frame containing two columns: `Z` is a treatment indicator, `B` is
a blocking factor
