# Create missingness indicator variables and non-informatively fill in missing values

Given a `data.frame` or `formula` and data, `fill.NAs()` returns an
expanded data frame, including a new missingness flag for each variable
with missing values and replacing each missing entry with a value
representing a reasonable default for missing values in its column.
Functions in the formula are supported, with transformations happening
before `NA` replacement. The expanded data frame is useful for
propensity modeling and balance checking when there are covariates with
missing values.

## Usage

``` r
fill.NAs(x, data = NULL, all.covs = FALSE, contrasts.arg = NULL)
```

## Arguments

- x:

  Can be either a data frame (in which case the data argument should be
  `NULL`) or a formula (in which case data must be a data.frame)

- data:

  If x is a formula, this must be a data.frame. Otherwise it will be
  ignored.

- all.covs:

  Should the response variable be imputed? For formula `x`, this is the
  variable on the left hand side. For `data.frame` `x`, the response is
  considered the first column.

- contrasts.arg:

  (from `model.matrix`) A list, whose entries are values (numeric
  matrices or character strings naming functions) to be used as
  replacement values for the
  [`contrasts`](https://rdrr.io/r/stats/contrasts.html) replacement
  function and whose names are the names of columns of `data` containing
  [`factor`](https://rdrr.io/r/base/factor.html)s.

## Value

A `data.frame` with all `NA` values replaced with mean values and
additional indicator columns for each column including missing values.
Suitable for directly passing to [`lm`](https://rdrr.io/r/stats/lm.html)
or other model building functions to build propensity scores.

## Details

`fill.NAs` prepares data for use in a model or matching procedure by
filling in missing values with minimally invasive substitutes. Fill-in
is performed column-wise, with each column being treated individually.
For each column that is missing, a new column is created of the form
“ColumnName.NA” with indicators for each observation that is missing a
value for “ColumnName”. Rosenbaum and Rubin (1984, Sec. 2.4 and Appendix
B) discuss propensity score models using this data structure.

The replacement value used to fill in a missing value is simple mean
replacement. For transformations of variables, e.g. `y ~ x1 * x2`, the
transformation occurs first. The transformation column will be `NA` if
any of the base columns are `NA`. Fill-in occurs next, replacing all
missing values with the observed column mean. This includes
transformation columns.

Data can be passed to `fill.NAs` in two ways. First, you can simply pass
a `data.frame` object and `fill.NAs` will fill every column.
Alternatively, you can pass a `formula` and a `data.frame`. Fill-in will
only be applied to columns specifically used in the formula. Prior to
fill-in, any functions in the formula will be expanded. If any arguments
to the functions are `NA`, the function value will also be `NA` and
subject to fill-in.

By default, `fill.NAs` does not impute the response variable. This is to
encourage more sophisticated imputation schemes when the response is a
treatment indicator in a matching problem. This behavior can be
overridden by setting `all.covs = TRUE`.

## References

Rosenbaum, Paul R. and Rubin, Donald B. (1984) ‘Reducing Bias in
Observational Studies using Subclassification on the Propensity Score,’
*Journal of the American Statistical Association*, **79**, 516 – 524.

Von Hipple, Paul T. (2009) ‘How to impute interactions, squares, and
other transformed variables,’ *Sociological Methodology*, **39**(1), 265
– 291.

## See also

[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`lm`](https://rdrr.io/r/stats/lm.html)

## Author

Mark M. Fredrickson and Jake Bowers

## Examples

``` r
data(nuclearplants)
### Extract some representative covariates:
np.missing <- nuclearplants[c('t1', 't2', 'ne', 'ct', 'cum.n')]

 ### create some missingness in the covariates
 n <- dim(np.missing)[1]
 k <- dim(np.missing)[2]

 for (i in 1:n) {
   missing <- rbinom(1, prob = .1, size = k)
   if (missing > 0) {
     np.missing[i, sample(k, missing)] <- NA
   }
 }

### Restore outcome and treatment variables:
np.missing <- data.frame(nuclearplants[c('cost', 'pr')], np.missing)

### Fit a propensity score but with missing covariate data flagged
### and filled in, as in Rosenbaum and Rubin (1984, Appendix):
np.filled <- fill.NAs(pr ~ t1 * t2, np.missing)
# Look at np.filled to establish what missingness flags were created
head(np.filled)
#>   pr       t1 t2    t1:t2 t1.NA t2.NA
#> H  0 13.72414 46 840.8333  TRUE FALSE
#> I  0 10.00000 73 730.0000 FALSE FALSE
#> A  1 10.00000 85 850.0000 FALSE FALSE
#> J  0 11.00000 67 737.0000 FALSE FALSE
#> B  1 11.00000 78 858.0000 FALSE FALSE
#> K  0 13.00000 51 663.0000 FALSE FALSE
(np.glm <- glm(pr ~ ., family=binomial, data=np.filled))
#> 
#> Call:  glm(formula = pr ~ ., family = binomial, data = np.filled)
#> 
#> Coefficients:
#> (Intercept)           t1           t2      `t1:t2`    t1.NATRUE    t2.NATRUE  
#>  -39.532987     0.594190     0.432671     0.002249    -0.583510     0.601206  
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  26 Residual
#> Null Deviance:       39.75 
#> Residual Deviance: 21.06     AIC: 33.06
(glm(pr ~ t1 + t2 + `t1:t2` + t1.NA + t2.NA,
                family=binomial, data=np.filled))
#> 
#> Call:  glm(formula = pr ~ t1 + t2 + `t1:t2` + t1.NA + t2.NA, family = binomial, 
#>     data = np.filled)
#> 
#> Coefficients:
#> (Intercept)           t1           t2      `t1:t2`    t1.NATRUE    t2.NATRUE  
#>  -39.532987     0.594190     0.432671     0.002249    -0.583510     0.601206  
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  26 Residual
#> Null Deviance:       39.75 
#> Residual Deviance: 21.06     AIC: 33.06
# In a non-interactive session, the following may help, as long as
# the formula passed to `fill.NAs` (plus any missingness flags) is
# the desired formula for the glm.
(glm(formula(terms(np.filled)), family=binomial, data=np.filled))
#> 
#> Call:  glm(formula = formula(terms(np.filled)), family = binomial, data = np.filled)
#> 
#> Coefficients:
#> (Intercept)           t1           t2      `t1:t2`    t1.NATRUE    t2.NATRUE  
#>  -39.532987     0.594190     0.432671     0.002249    -0.583510     0.601206  
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  26 Residual
#> Null Deviance:       39.75 
#> Residual Deviance: 21.06     AIC: 33.06

### produce a matrix of propensity distances based on the propensity model
### with fill-in and flagging. Then perform pair matching on it:
pairmatch(match_on(np.glm, data=np.filled), data=np.filled)
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.4  1.1  1.3  1.2 <NA> <NA> <NA>  1.3  1.8 <NA> <NA>  1.6  1.2 1.10 <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.9  1.4  1.7  1.5 <NA>  1.6  1.5  1.7 <NA>  1.1 <NA> <NA> <NA>  1.8  1.9 1.10 

## fill NAs without using treatment contrasts by making a list of contrasts for
## each factor ## following hints from https://stackoverflow.com/a/4569239/161808

np.missing$t1F<-factor(np.missing$t1)
cov.factors <- sapply(np.missing[,c("t1F","t2")],is.factor)
cov.contrasts <- lapply(
  np.missing[,names(cov.factors)[cov.factors],drop=FALSE],
  contrasts, contrasts = FALSE)

## make a data frame filling the missing covariate values, but without
## excluding any levels of any factors
np.noNA2<-fill.NAs(pr~t1F+t2,data=np.missing,contrasts.arg=cov.contrasts)
```
