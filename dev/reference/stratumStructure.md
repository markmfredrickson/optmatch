# Return structure of matched sets

Tabulate treatment:control ratios occurring in matched sets, and the
frequency of their occurrence.

## Usage

``` r
stratumStructure(stratum, trtgrp = NULL, min.controls = 0, max.controls = Inf)

# S3 method for class 'optmatch'
stratumStructure(stratum, trtgrp, min.controls = 0, max.controls = Inf)

# Default S3 method
stratumStructure(stratum, trtgrp, min.controls = 0, max.controls = Inf)

# S3 method for class 'stratumStructure'
print(x, ...)
```

## Arguments

- stratum:

  Matched strata, as returned by
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

- trtgrp:

  Dummy variable for treatment group membership. (Not required if
  `stratum` is an optmatch object, as returned by
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).)

- min.controls:

  For display, the number of treatment group members per stratum will be
  truncated at the reciprocal of `min.controls`.

- max.controls:

  For display, the number of control group members will be truncated at
  `max.controls`.

- x:

  stratumStructure object to be printed.

- ...:

  Additional arguments to `print`.

## Value

A table showing frequency of occurrence of those treatment:control
ratios that occur.

The ‘effective sample size’ of the stratification, in matched pairs.
Given as an attribute of the table, named
‘`comparable.num.matched.pairs`’; see Note.

## Note

For comparing treatment and control groups both of size 10, say, a
stratification consisting of two strata, one with 9 treatments and 1
control, has a smaller ‘effective sample size’, intuitively, than a
stratification into 10 matched pairs, despite the fact that both contain
20 subjects in total. `stratumStructure` first summarizes this aspect of
the structure of the stratification it is given, then goes on to
identify one number as the stratification's effective sample size. The
‘`comparable.num.matched.pairs`’ attribute returned by
`stratumStructure` is the sum of harmonic means of the sizes of the
treatment and control subgroups of each stratum, a general way of
calibrating such differences as well as differences in the number of
subjects contained in a stratification. For example, by this metric the
9:1, 1:9 stratification is comparable to 3.6 matched pairs.

Why should effective sample size be calculated this way? The phrase
‘effective sample size’ suggests the observations are taken to be
similar in information content. Modeling them as random variables, this
suggests that they be assumed to have the same variance, \\\sigma\\,
conditional on what stratum they reside in. If that is the case, and if
also treatment and control observations differ in expectation by a
constant that is the same for each stratum, then it can be shown that
the optimum weights with which to combine treatment-control contrasts
across strata, \\s\\, are proportional to the stratum-wise harmonic
means of treatment and control counts, \\h_s = \[(n\_{ts}^{-1} +
n\_{cs}^{-1})/2\]^{-1}\\ (Kalton, 1968). The thus-weighted average of
contrasts then has variance \\2\sigma/\sum_s h_s\\. This motivates the
use of \\\sum_s h_s\\ as a measure of effective sample size (Hansen,
2011). Somewhat different motivations of the same calculation appear in
Hansen (2004) and in Hansen and Bowers (2008). Since for a matched pair
\\s\\, \\h_s=1\\, \\\sum_s h_s\\ can be thought of as the number of
matched pairs needed to attain comparable precision.

## References

Kalton, G. (1968), ‘Standardization: A technique to control for
extraneous variables’, *Applied Statistics*, **17**, 118–136.

Hansen, B.B. (2004), ‘Full Matching in an Observational Study of
Coaching for the SAT’, *Journal of the American Statistical
Association*, **99**, 609–618.

Hansen B.B. and Bowers, J. (2008), ‘Covariate balance in simple,
stratified and clustered comparative studies’, *Statistical Science*,
**23** (2), 219–236.

Hansen, B.B. (2011), ‘Propensity score matching to extract latent
experiments from nonexperimental data: A case study’. Ch. 9 of *Looking
Backwards: Proceedings from a Conference in Honor of Paul W. Holland*,
Springer.

## See also

[`matched`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)

## Author

Ben B. Hansen

## Examples

``` r
data(plantdist)
plantsfm <- fullmatch(plantdist) # A full match with unrestricted
#> Warning: Without 'data' argument the order of the match is not guaranteed
#>     to be the same as your original data.
                                 # treatment-control balance
plantsfm1 <- fullmatch(plantdist,min.controls=2, max.controls=3)
#> Warning: Without 'data' argument the order of the match is not guaranteed
#>     to be the same as your original data.
stratumStructure(plantsfm)
#> 1:1 1:2 1:3 1:4 1:6 
#>   2   2   1   1   1 
stratumStructure(plantsfm1)
#> 1:2 1:3 
#>   2   5 
stratumStructure(plantsfm, max.controls=4)
#>  1:1  1:2  1:3 1:4+ 
#>    2    2    1    2 
```
