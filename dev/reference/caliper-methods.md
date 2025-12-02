# Prepare matching distances suitable for matching within calipers.

Encodes calipers, or maximum allowable distances within which to match.
The result of a call to `caliper` is itself a distance specification
between treated and control units that can be used with
[`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
or
[`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md).
Calipers can also be combined with other distance specifications for
richer matching problems.

## Usage

``` r
caliper(x, width, exclude = c(), compare = `<=`, values = FALSE)

# S4 method for class 'InfinitySparseMatrix'
caliper(x, width, exclude = c(), compare = `<=`, values = FALSE)

# S4 method for class 'matrix'
caliper(x, width, exclude = c(), compare = `<=`, values = FALSE)

# S4 method for class 'optmatch.dlist'
caliper(x, width, exclude = c(), compare = `<=`, values = FALSE)
```

## Arguments

- x:

  A distance specification created with
  [`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  or similar.

- width:

  The width of the caliper: how wide of a margin to allow in matches. Be
  careful in setting the width. Vector valued arguments will be recycled
  for each of the finite entries in `x` (and no order is guaranteed for
  `x` for some types of distance objects).

- exclude:

  (Optional) A character vector of observations (corresponding to row
  and column names) to exclude from the caliper.

- compare:

  A function that decides that whether two observations are with the
  caliper. The default is `` `<=` ``. `` `<` `` is a common alternative.

- values:

  Should the returned object be made of all zeros (`values = FALSE`, the
  default) or should the object include the values of the original
  object (`values = TRUE`)?

## Value

A matrix like object that is suitable to be given as `distance` argument
to
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
or
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).
The caliper will be only zeros and `Inf` values, indicating a possible
match or no possible match, respectively.

You can combine the results of `caliper` with other distances using the
`` `+` `` operator. See the examples for usage.

## Details

`caliper` is a generic function with methods for any of the allowed
distance specifications: user created matrices, the results of
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
the results of
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
or combinations (using `` `+` ``) of these objects.

`width` provides the size of the caliper, the allowable distance for
matching. If the distance between a treated and control pair is less
than or equal to this distance, it is allowed kept; otherwise, the pair
is discarded from future matching. The default comparison of "equal or
less than can" be changed to any other comparison function using the
`comparison` argument.

It is important to understand that `width` argument is defined on the
scale of these distances. For univariate distances such as propensity
scores, it is common to specify calipers in standard deviations. If a
caliper of this nature is desired, you must either find the standard
deviation directly or use the
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
function with its `caliper` argument. Since `match_on` has access to the
underlying univariate scores, for example for the GLM method, it can
determine the caliper width in standard deviations.

If you wish to exclude specific units from the caliper requirements,
pass the names of these units in the `exclude` argument. These units
will be allowed to match any other unit.

## References

P.~R. Rosenbaum and D.~B. Rubin (1985), ‘Constructing a control group
using multivariate matched sampling methods that incorporate the
propensity score’, *The American Statistician*, **39** 33–38.

## See also

[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

## Author

Mark M. Fredrickson and Ben B. Hansen

## Examples

``` r
data(nuclearplants)


### Caliper of 100 MWe on plant capacity
caliper(match_on(pr~cap, data=nuclearplants, method="euclidean"), width=100)
#>          control
#> treatment   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X
#>         A Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         B Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         C Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         D Inf Inf Inf   0 Inf   0 Inf   0 Inf   0 Inf Inf Inf Inf Inf Inf Inf
#>         E Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         F Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         G Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         a Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0 Inf   0 Inf   0   0
#>         b   0 Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf Inf   0   0 Inf Inf   0
#>         c Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0 Inf   0 Inf   0   0
#>          control
#> treatment   Y   Z   d   e   f
#>         A Inf   0 Inf Inf Inf
#>         B Inf   0 Inf Inf Inf
#>         C Inf Inf   0   0   0
#>         D   0 Inf Inf Inf Inf
#>         E Inf   0 Inf Inf Inf
#>         F Inf Inf   0   0   0
#>         G Inf Inf   0   0   0
#>         a Inf Inf Inf   0   0
#>         b Inf Inf   0   0 Inf
#>         c Inf Inf Inf   0   0

### Caliper of 1/2 a pooled SD of plant capacity
caliper(match_on(pr~cap, data=nuclearplants), width=.5)
#>          control
#> treatment   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X
#>         A Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         B Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         C Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         D Inf Inf Inf   0 Inf   0 Inf   0 Inf   0 Inf Inf Inf Inf Inf Inf Inf
#>         E Inf   0   0 Inf Inf Inf Inf Inf Inf Inf   0 Inf Inf Inf   0 Inf Inf
#>         F Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         G Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf   0   0   0 Inf   0   0
#>         a Inf Inf Inf Inf   0 Inf   0 Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf
#>         b   0 Inf Inf Inf   0 Inf   0 Inf   0 Inf Inf Inf   0 Inf Inf Inf   0
#>         c Inf Inf Inf Inf   0 Inf   0 Inf Inf Inf Inf   0 Inf   0 Inf   0 Inf
#>          control
#> treatment   Y   Z   d   e   f
#>         A Inf   0 Inf Inf Inf
#>         B Inf   0 Inf Inf Inf
#>         C Inf Inf   0   0   0
#>         D   0 Inf Inf Inf Inf
#>         E Inf   0 Inf Inf Inf
#>         F Inf Inf   0   0   0
#>         G Inf Inf   0   0   0
#>         a Inf Inf Inf   0   0
#>         b Inf Inf   0   0 Inf
#>         c Inf Inf Inf   0   0

### Caliper  of .2 pooled SDs in the propensity score
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
ppty.dist <- match_on(ppty)

pptycaliper <- caliper(ppty.dist, width = .2)

### caliper on the Mahalanobis distance
caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3)
#>          control
#> treatment   H   I   J   K   L M N   O   P Q   R S   T   U V   W X   Y   Z d   e
#>         A Inf   0   0 Inf Inf 0 0   0 Inf 0 Inf 0   0   0 0   0 0 Inf Inf 0 Inf
#>         B Inf   0   0 Inf Inf 0 0   0 Inf 0   0 0   0   0 0   0 0 Inf   0 0 Inf
#>         C   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0   0
#>         D Inf   0   0 Inf Inf 0 0   0 Inf 0   0 0   0   0 0   0 0 Inf   0 0 Inf
#>         E   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         F Inf Inf Inf Inf Inf 0 0 Inf Inf 0   0 0 Inf   0 0 Inf 0   0   0 0 Inf
#>         G   0   0   0   0 Inf 0 0   0 Inf 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         a   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         b Inf   0   0   0 Inf 0 0   0 Inf 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         c Inf   0   0 Inf Inf 0 0   0 Inf 0 Inf 0   0 Inf 0   0 0 Inf Inf 0 Inf
#>          control
#> treatment   f
#>         A   0
#>         B   0
#>         C   0
#>         D   0
#>         E   0
#>         F Inf
#>         G   0
#>         a   0
#>         b   0
#>         c   0

### Combining a Mahalanobis distance matching with a caliper
### of 1 pooled SD in the propensity score:
mhd.pptyc <- caliper(ppty.dist, width = 1) +
          match_on(pr ~ t1 + t2, data = nuclearplants)
pairmatch(mhd.pptyc, data = nuclearplants)
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> 

### Excluding observations from caliper requirements:
caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3, exclude = c("A", "f"))
#>          control
#> treatment   H   I   J   K   L M N   O   P Q   R S   T   U V   W X   Y   Z d   e
#>         A   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0   0
#>         B Inf   0   0 Inf Inf 0 0   0 Inf 0   0 0   0   0 0   0 0 Inf   0 0 Inf
#>         C   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0   0
#>         D Inf   0   0 Inf Inf 0 0   0 Inf 0   0 0   0   0 0   0 0 Inf   0 0 Inf
#>         E   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         F Inf Inf Inf Inf Inf 0 0 Inf Inf 0   0 0 Inf   0 0 Inf 0   0   0 0 Inf
#>         G   0   0   0   0 Inf 0 0   0 Inf 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         a   0   0   0   0   0 0 0   0   0 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         b Inf   0   0   0 Inf 0 0   0 Inf 0   0 0   0   0 0   0 0   0   0 0 Inf
#>         c Inf   0   0 Inf Inf 0 0   0 Inf 0 Inf 0   0 Inf 0   0 0 Inf Inf 0 Inf
#>          control
#> treatment f
#>         A 0
#>         B 0
#>         C 0
#>         D 0
#>         E 0
#>         F 0
#>         G 0
#>         a 0
#>         b 0
#>         c 0

### Returning values directly (equal up to the the attributes)
all(abs((caliper(ppty.dist, 1) + ppty.dist) -
        caliper(ppty.dist, 1, values = TRUE)) < .Machine$Double.eps)
#> [1] TRUE
```
