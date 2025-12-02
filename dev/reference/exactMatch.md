# Generate an exact matching set of subproblems.

An exact match is one based on a factor. Within a level, all
observations are allowed to be matched. An exact match can be combined
with another distance matrix to create a set of matching subproblems.

## Usage

``` r
exactMatch(x, ...)

# S4 method for class 'vector'
exactMatch(x, treatment)

# S4 method for class 'formula'
exactMatch(x, data = NULL, subset = NULL, na.action = NULL, ...)
```

## Arguments

- x:

  A factor vector or a formula, used to select method.

- ...:

  Additional arguments for methods.

- treatment:

  A logical or binary vector the same length as `x` indicating treatment
  and control for each unit in the study. TRUE or 1 represents a
  treatment unit, FALSE or 0 represents a control unit. NA units are
  excluded.

- data:

  A `data.frame` or `matrix` that contains the variables used in the
  formula `x`.

- subset:

  an optional vector specifying a subset of observations to be used

- na.action:

  A function which indicates what should happen when the data contain
  `NA`s

## Value

A matrix like object, which is suitable to be given as `distance`
argument to
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
or
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).
The exact match will be only zeros and `Inf` values, indicating a
possible match or no possible match, respectively. It can be added to a
another distance matrix to create a subclassed matching problem.

## Details

`exactMatch` creates a block diagonal matrix of 0s and `Inf`s. The pairs
with 0 entries are within the same level of the factor and legitimate
matches. `Inf` indicates units in different levels. `exactMatch`
replaces the `structure.fmla` argument to several functions in previous
versions of optmatch. For the `factor` method, the two vectors `x` and
`treatment` must be the same length. The vector `x` is interpreted as
indicating the grouping factors for the data, and the vector `treatment`
indicates whether a unit is in the treatment or control groups. At least
one of these two vectors must have names. For the `formula` method, the
`data` argument may be omitted, in which case the method attempts to
find the variables in the environment from which the function was
called. This behavior, and the arguments `subset` and `na.action`,
mimics the behavior of [`lm`](https://rdrr.io/r/stats/lm.html).

## See also

[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
[`antiExactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/antiExactMatch.md),
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

## Author

Mark M. Fredrickson

## Examples

``` r
data(nuclearplants)

### First generate a standard propensity score
ppty <- glm(pr~.-(pr+cost), family = binomial(), data = nuclearplants)
ppty.distances <- match_on(ppty)

### Only allow matches within the partial turn key plants
pt.em <- exactMatch(pr ~ pt, data = nuclearplants)
as.matrix(pt.em)
#>          control
#> treatment   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X
#>         A   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         B   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         C   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         D   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         E   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         F   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         G   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#>         a Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
#>         b Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
#>         c Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
#>          control
#> treatment   Y   Z   d   e   f
#>         A   0   0 Inf Inf Inf
#>         B   0   0 Inf Inf Inf
#>         C   0   0 Inf Inf Inf
#>         D   0   0 Inf Inf Inf
#>         E   0   0 Inf Inf Inf
#>         F   0   0 Inf Inf Inf
#>         G   0   0 Inf Inf Inf
#>         a Inf Inf   0   0   0
#>         b Inf Inf   0   0   0
#>         c Inf Inf   0   0   0

### Blunt matches:
match.pt.em <- fullmatch(pt.em)
#> Warning: Without 'data' argument the order of the match is not guaranteed
#>     to be the same as your original data.
print(match.pt.em, grouped = TRUE)
#>  Group                                  Members
#>    0.1                                     A, Z
#>    0.2                                     B, Y
#>    0.3                                     C, X
#>    0.4                                     D, W
#>    0.5                                     E, V
#>    0.6                                     F, U
#>    0.7 G, H, I, J, K, L, M, N, O, P, Q, R, S, T
#>    1.1                                     a, f
#>    1.2                                     b, e
#>    1.3                                     c, d

### Combine the propensity scores with the subclasses:
match.ppty.em <- fullmatch(ppty.distances + pt.em)
#> Warning: Without 'data' argument the order of the match is not guaranteed
#>     to be the same as your original data.
print(match.ppty.em, grouped = TRUE)
#>  Group                                        Members
#>    0.1                                     A, D, F, V
#>    0.2                                        B, G, U
#>    0.3 C, H, J, K, L, M, N, O, P, Q, R, T, W, X, Y, Z
#>    0.5                                        E, I, S
#>    1.1                                           a, d
#>    1.2                                           b, f
#>    1.3                                           c, e
```
