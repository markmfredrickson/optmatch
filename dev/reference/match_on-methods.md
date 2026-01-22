# Create treated to control distances for matching problems

A function with which to produce matching distances, for instance
Mahalanobis distances, propensity score discrepancies or calipers, or
combinations thereof, for
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
or
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
to subsequently “match on”. Conceptually, the result of a call
`match_on` is a treatment-by-control matrix of distances. Because these
matrices can grow quite large, in practice `match_on` produces either an
ordinary dense matrix or a special sparse matrix structure (that can
make use of caliper and exact matching constraints to reduce storage
requirements). Methods are supplied for these sparse structures,
`InfinitySparseMatrix`es, so that they can be manipulated and modified
in much the same way as dense matrices.

## Usage

``` r
match_on(x, within = NULL, caliper = NULL, exclude = NULL, data = NULL, ...)

# S3 method for class 'glm'
match_on(
  x,
  within = NULL,
  caliper = NULL,
  exclude = NULL,
  data = NULL,
  standardization.scale = NULL,
  ...
)

# S3 method for class 'bigglm'
match_on(
  x,
  within = NULL,
  caliper = NULL,
  exclude = NULL,
  data = NULL,
  standardization.scale = NULL,
  ...
)

# S3 method for class 'formula'
match_on(
  x,
  within = NULL,
  caliper = NULL,
  exclude = NULL,
  data = NULL,
  subset = NULL,
  method = "mahalanobis",
  ...
)

# S3 method for class '`function`'
match_on(
  x,
  within = NULL,
  caliper = NULL,
  exclude = NULL,
  data = NULL,
  z = NULL,
  ...
)

# S3 method for class 'numeric'
match_on(x, within = NULL, caliper = NULL, exclude = NULL, data = NULL, z, ...)

# S3 method for class 'InfinitySparseMatrix'
match_on(x, within = NULL, caliper = NULL, exclude = NULL, data = NULL, ...)

# S3 method for class 'matrix'
match_on(x, within = NULL, caliper = NULL, exclude = NULL, data = NULL, ...)
```

## Arguments

- x:

  A model formula, fitted glm or other object implicitly specifying a
  distance; see blurbs on specific methods in Details.

- within:

  A valid distance specification, such as the result of
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
  or
  [`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md).
  Finite entries indicate which distances to create. Including this
  argument can significantly speed up computation for sparse matching
  problems. Specify this filter either via `within` or via `strata`
  elements of a formula; mixing these methods will fail.

- caliper:

  The width of a caliper to use to exclude treated-control pairs with
  values greater than the width. For some methods, there may be a speed
  advantage to passing a width rather than using the
  [`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
  function on an existing distance specification.

- exclude:

  A list of units (treated or control) to exclude from the `caliper`
  argument (if supplied).

- data:

  An optional data frame.

- ...:

  Other arguments for methods.

- standardization.scale:

  Function for rescaling of `scores(x)`, or `NULL`; defaults to `mad`.
  (See Details.)

- subset:

  A subset of the data to use in creating the distance specification.

- method:

  A string indicating which method to use in computing the distances
  from the data. The current possibilities are
  `"mahalanobis", "euclidean"` or `"rank_mahalanobis"`.

- z:

  A logical or binary vector indicating treatment and control for each
  unit in the study. TRUE or 1 represents a treatment unit, FALSE of 0
  represents a control unit. Any unit with NA treatment status will be
  excluded from the distance matrix.

## Value

A distance specification (a matrix or similar object) which is suitable
to be given as the `distance` argument to
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
or
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

## Details

`match_on` is generic. There are several supplied methods, all providing
the same basic output: a matrix (or similar) object with treated units
on the rows and control units on the columns. Each cell \[i,j\] then
indicates the distance from a treated unit i to control unit j. Entries
that are `Inf` are said to be unmatchable. Such units are guaranteed to
never be in a matched set. For problems with many `Inf` entries, so
called sparse matching problems, `match_on` uses a special data type
that is more space efficient than a standard R `matrix`. When problems
are not sparse (i.e. dense), `match_on` uses the standard `matrix` type.

`match_on` methods differ on the types of arguments they take, making
the function a one-stop location of many different ways of specifying
matches: using functions, formulas, models, and even simple scores. Many
of the methods require additional arguments, detailed below. All methods
take a `within` argument, a distance specification made using
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
or
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
(or some additive combination of these or other distance creating
functions). All `match_on` methods will use the finite entries in the
`within` argument as a guide for producing the new distance. Any entry
that is `Inf` in `within` will be `Inf` in the distance matrix returned
by `match_on`. This argument can reduce the processing time needed to
compute sparse distance matrices.

Details for each particular first type of argument follow:

**First argument (`x`): `glm`.** The model is assumed to be a fitted
propensity score model. From this it extracts distances on the *linear*
propensity score: fitted values of the linear predictor, the link
function applied to the estimated conditional probabilities, as opposed
to the estimated conditional probabilities themselves (Rosenbaum &
Rubin, 1985). For example, a logistic model (`glm` with
`family=binomial()`) has the logit function as its link, so from such
models `match_on` computes distances in terms of logits of the estimated
conditional probabilities, i.e. the estimated log odds.

Optionally these distances are also rescaled. The default is to rescale,
by the reciprocal of an outlier-resistant variant of the pooled s.d. of
propensity scores; see
[`standardization_scale`](https://markmfredrickson.github.io/optmatch/dev/reference/standardization_scale.md).
(The `standardization.scale` argument of this function can be used to
change how this dispersion is calculated, e.g. to calculate an ordinary
not an outlier-resistant s.d.; it will be passed down to
`standardization_scale` as its `standardizer` argument.) To skip
rescaling, set argument `standardization.scale` to 1. The overall result
records absolute differences between treated and control units on
linear, possibly rescaled, propensity scores.

In addition, one can impose a caliper in terms of these distances by
providing a scalar as a `caliper` argument, forbidding matches between
treatment and control units differing in the calculated propensity score
by more than the specified caliper. For example, Rosenbaum and Rubin's
(1985) caliper of one-fifth of a pooled propensity score s.d. would be
imposed by specifying `caliper=.2`, in tandem either with the default
rescaling or, to follow their example even more closely, with the
additional specification `standardization.scale=sd`. Propensity calipers
are beneficial computationally as well as statistically, for reasons
indicated in the below discussion of the `numeric` method.

One can also specify exactMatching criteria by using `strata(foo)`
inside the formula to build the `glm`. For example, passing
`glm(y ~ x + strata(s))` to `match_on` is equivalent to passing
`within=exactMatch(y ~ strata(s))`. Note that when combining with the
`caliper` argument, the standard deviation used for the caliper will be
computed across all strata, not within each strata.

If data used to fit the glm have missing values in the left-hand side
(dependent) variable, these observations are omitted from the output of
match_on. If there are observations with missing values in right hand
side (independent) variables, then a re-fit of the model after imputing
these variables using a simple scheme and adding indicator variables of
missingness will be attempted, via the
[`scores`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
function.

**First argument (`x`): `bigglm`.** This method works analogously to the
`glm` method, but with `bigglm` objects, created by the `bigglm`
function from package ‘biglm’, which can handle bigger data sets than
the ordinary glm function can.

**First argument (`x`): `formula`.** The formula must have `Z`, the
treatment indicator (`Z=0` indicates control group, `Z=1` indicates
treatment group), on the left hand side, and any variables to compute a
distance on on the right hand side. E.g. `Z ~ X1 + X2`. The Mahalanobis
distance is calculated as the square root of d'Cd, where d is the vector
of X-differences on a pair of observations and C is an inverse
(generalized inverse) of the pooled covariance of Xes. (The pooling is
of the covariance of X within the subset defined by `Z==0` and within
the complement of that subset. This is similar to a Euclidean distance
calculated after reexpressing the Xes in standard units, such that the
reexpressed variables all have pooled SDs of 1; except that it addresses
redundancies among the variables by scaling down variables contributions
in proportion to their correlations with other included variables.)

Euclidean distance is also available, via `method="euclidean"`, as are
two flavors of ranked-based Mahalanobis distance, via
`method="rank_mahalanobis"` or `method="pooled_cov_rank_mahalanobis"`.
Either rank-transforms the covariates first; they differ in whether
subsequent covariance of thus-transformed covariates is calculated on
all subjects or by pooling of with-group covariances across treatment
and control. The `method=` argument can be abbreviated in the usual way
(via \[base::pmatch()\]).

The treatment indicator `Z` as noted above must either be numeric (1
representing treated units and 0 control units) or logical (`TRUE` for
treated, `FALSE` for controls). (Earlier versions of the software
accepted factor variables and other types of numeric variable; you may
have to update existing scripts to get them to run.)

As an alternative to specifying a `within` argument, when `x` is a
formula, the `strata` command can be used inside the formula to specify
exact matching. For example, rather than using
`within=exactMatch(y ~ z, data=data)`, you may update your formula as
`y ~ x + strata(z)`. Do not use both methods (`within` and `strata`
simultaneously. Note that when combining with the `caliper` argument,
the standard deviation used for the caliper will be computed across all
strata, not separately by stratum.

A unit with NA treatment status (`Z`) is ignored and will not be
included in the distance output. Missing values in variables on the
right hand side of the formula are handled as follows. By default
`match_on` will (1) create a matrix of distances between observations
which have only valid values for \*\*all\*\* covariates and then (2)
append matrices of Inf values for distances between observations either
of which has a missing values on any of the right-hand-side variables.
(I.e., observations with missing values are retained in the output, but
matches involving them are forbidden.)

**First argument (`x`): `function`.** The passed function must take
arguments: `index`, `data`, and `z`. The `data` and `z` arguments will
be the same as those passed directly to `match_on`. The `index` argument
is a matrix of two columns, representing the pairs of treated and
control units that are valid comparisons (given any `within` arguments).
The first column is the row name or id of the treated unit in the `data`
object. The second column is the id for the control unit, again in the
`data` object. For each of these pairs, the function should return the
distance between the treated unit and control unit. This may sound
complicated, but is simple to use. For example, a function that returned
the absolute difference between two units using a vector of data would
be
` f <- function(index, data, z) { abs(data[index[,1]] - data[index[,2]]) } `.
(Note: This simple case is precisely handled by the `numeric` method.)

**First argument (`x`): `numeric`.** This returns absolute differences
between treated and control units' values of `x`. If a caliper is
specified, pairings with `x`-differences greater than it are forbidden.
Conceptually, those distances are set to `Inf`; computationally, if
either of `caliper` and `within` has been specified then only
information about permissible pairings will be stored, so the forbidden
pairings are simply omitted. Providing a `caliper` argument here, as
opposed to omitting it and afterward applying the
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
function, reduces storage requirements and may otherwise improve
performance, particularly in larger problems.

For the numeric method, `x` must have names. If `z` is named it must
have the same names as `x`, though it allows for a different ordering of
names. `x`'s name ordering is considered canonical.

**First argument (`x`): `matrix` or `InfinitySparseMatrix`.** These just
return their arguments as these objects are already valid distance
specifications.

## References

P.~R. Rosenbaum and D.~B. Rubin (1985), ‘Constructing a control group
using multivariate matched sampling methods that incorporate the
propensity score’, *The American Statistician*, **39** 33–38.

## See also

[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md),
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)

[`scores`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)

## Examples

``` r
data(nuclearplants)
match_on.examples <- list()
### Propensity score distances.
### Recommended approach:
(aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants))
#> 
#> Call:  glm(formula = pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
#> 
#> Coefficients:
#> (Intercept)         date           t1           t2          cap           ne  
#>   43.909419    -1.220088     0.938645     0.396002     0.001197    -0.995651  
#>          ct           bw        cum.n           pt  
#>   -2.615671    -0.088696     0.033575    -0.352112  
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  22 Residual
#> Null Deviance:       39.75 
#> Residual Deviance: 20.92     AIC: 40.92
match_on.examples$ps1 <- match_on(aGlm)
### A second approach: first extract propensity scores, then separately
### create a distance from them.  (Useful when importing propensity
### scores from an external program.)
plantsPS <- predict(aGlm)
match_on.examples$ps2 <- match_on(pr~plantsPS, data=nuclearplants)
### Full matching on the propensity score.
fm1 <- fullmatch(match_on.examples$ps1, data = nuclearplants)
fm2 <- fullmatch(match_on.examples$ps2, data = nuclearplants)
### Because match_on.glm uses robust estimates of spread,
### the results differ in detail -- but they are close enough
### to yield similar optimal matches.
all(fm1 == fm2) # The same
#> [1] TRUE

### Mahalanobis distance:
match_on.examples$mh1 <- match_on(pr ~ t1 + t2, data = nuclearplants)

### Absolute differences on a scalar:
tmp <- nuclearplants$t1
names(tmp) <- rownames(nuclearplants)

(absdist <- match_on(tmp, z = nuclearplants$pr,
                  within = exactMatch(pr ~ pt, nuclearplants)))
#> $`0`
#>          control
#> treatment H  I  J K  L M  N  O P Q R S  T U V  W X Y  Z
#>         A 4  0  1 3  2 4  2  2 3 5 7 3  1 8 5  1 6 9 10
#>         B 3  1  0 2  1 3  1  1 2 4 6 2  0 7 4  0 5 8  9
#>         C 1  5  4 2  3 1  3  3 2 0 2 2  4 3 0  4 1 4  5
#>         D 1  5  4 2  3 1  3  3 2 0 2 2  4 3 0  4 1 4  5
#>         E 2  6  5 3  4 2  4  4 3 1 1 3  5 2 1  5 0 3  4
#>         F 8 12 11 9 10 8 10 10 9 7 5 9 11 4 7 11 6 3  2
#>         G 5  9  8 6  7 5  7  7 6 4 2 6  8 1 4  8 3 0  1
#> 
#> $`1`
#>          control
#> treatment d e f
#>         a 1 3 0
#>         b 0 4 1
#>         c 6 2 5
#> 

### Pair matching on the variable `t1`:
pairmatch(absdist, data = nuclearplants)
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  0.1  0.1 <NA>  0.2 <NA> <NA> <NA>  0.3 <NA> <NA> <NA>  0.4 <NA> <NA> <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#> <NA>  0.4  0.3  0.5  0.2  0.6  0.5  0.7  0.7  0.6  1.2  1.3  1.1  1.1  1.2  1.3 


### Propensity score matching within subgroups:
match_on.examples$ps3 <- match_on(aGlm, exactMatch(pr ~ pt, nuclearplants))
fullmatch(match_on.examples$ps3, data = nuclearplants)
#>   H   I   A   J   B   K   L   M   C   N   O   P   Q   R   S   T   U   D   V   E 
#> 0.3 0.5 0.1 0.3 0.2 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.5 0.3 0.2 0.1 0.1 0.5 
#>   W   F   X   G   Y   Z   d   e   f   a   b   c 
#> 0.3 0.1 0.3 0.2 0.3 0.3 1.1 1.3 1.2 1.1 1.2 1.3 

### Propensity score matching with a propensity score caliper:
match_on.examples$pscal <- match_on.examples$ps1 + caliper(match_on.examples$ps1, 1)
fullmatch(match_on.examples$pscal, data = nuclearplants) # Note that the caliper excludes some units
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.4 <NA>  1.2  1.1 <NA> <NA>  1.2  1.2  1.2  1.2 <NA>  1.2  1.2  1.4  1.2 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.1  1.3  1.3  1.4  1.2  1.3  1.2  1.1 <NA>  1.2  1.4 <NA>  1.2  1.3  1.3  1.1 

### A Mahalanobis distance for matching within subgroups:
match_on.examples$mh2 <- match_on(pr ~ t1 + t2 , data = nuclearplants,
                            within = exactMatch(pr ~ pt, nuclearplants))

### Mahalanobis matching within subgroups, with a propensity score
### caliper:
fullmatch(match_on.examples$mh2 + caliper(match_on.examples$ps3, 1), data = nuclearplants)
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  0.1 <NA>  0.2  0.1 <NA> <NA>  0.2  0.2  0.2  0.2 <NA>  0.2  0.2  0.4  0.2 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  0.5  0.3  0.3  0.4  0.2  0.5  0.2  0.5 <NA>  0.4  1.1 <NA> <NA>  1.1 <NA>  1.1 

### Alternative methods to matching without groups (exact matching)
m1 <- match_on(pr ~ t1 + t2, data=nuclearplants, within=exactMatch(pr ~ pt, nuclearplants))
m2 <- match_on(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
# m1 and m2 are identical

m3 <- match_on(glm(pr ~ t1 + t2 + cost, data=nuclearplants,
                   family=binomial),
               data=nuclearplants,
               within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- match_on(glm(pr ~ t1 + t2 + cost + pt, data=nuclearplants,
                   family=binomial),
               data=nuclearplants,
               within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- match_on(glm(pr ~ t1 + t2 + cost + strata(pt), data=nuclearplants,
                   family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
```
