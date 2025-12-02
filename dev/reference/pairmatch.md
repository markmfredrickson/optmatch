# Optimal 1:1 and 1:k matching

Given a treatment group, a larger control reservoir, and a method for
creating discrepancies between each treatment and control unit (or
optionally an already created such discrepancy matrix), finds a pairing
of treatment units to controls that minimizes the sum of discrepancies.

## Usage

``` r
pairmatch(x, controls = 1, data = NULL, remove.unmatchables = FALSE, ...)

pair(x, controls = 1, data = NULL, remove.unmatchables = FALSE, ...)
```

## Arguments

- x:

  Any valid input to `match_on`. If `x` is a numeric vector, there must
  also be passed a vector `z` indicating grouping. Both vectors must be
  named.

  Alternatively, a precomputed distance may be entered.

- controls:

  The number of controls to be matched to each treatment

- data:

  Optional data set.

- remove.unmatchables:

  Should treatment group members for which there are no eligible
  controls be removed prior to matching?

- ...:

  Additional arguments to pass to
  [`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  (e.g. `within`)) or to
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  (e.g. `tol`). It is an error to pass `min.controls`, `max.controls`,
  `mean.controls` or `omit.fraction` as `pairmatch` must set these
  values.

## Value

A
[`optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
object (`factor`) indicating matched groups.

## Details

This is a wrapper to
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md);
see its documentation for more information, especially on additional
arguments to pass, additional discussion of valid input for parameter
`x`, and feasibility recovery.

If `remove.unmatchables` is `FALSE`, then if there are unmatchable
treated units then the matching as a whole will fail and no units will
be matched. If `TRUE`, then this unit will be removed and the function
will attempt to match each of the other treatment units. As of version
0.9-8, if there are fewer matchable treated units than matchable
controls then `pairmatch` will attempt to place each into a matched pair
each of the matchable controls and a strict subset of the matchable
treated units. (Previously matching would have failed for subclasses of
this structure.)

Matching can still fail, even with `remove.unmatchables` set to `TRUE`,
if there is too much competition for certain controls; if you find
yourself in that situation you should consider full matching, which
necessarily finds a match for everyone with an eligible match somewhere.

The units of the `optmatch` object returned correspond to members of the
treatment and control groups in reference to which the matching problem
was posed, and are named accordingly; the names are taken from the row
and column names of `distance` (with possible additions from the
optional `data` argument). Each element of the vector is the
concatenation of: (i) a character abbreviation of `subclass.indices`, if
that argument was given, or the string '`m`' if it was not; (ii) the
string `.`; and (iii) a non-negative integer. Unmatched units have `NA`
entries. Secondarily, `fullmatch` returns various data about the
matching process and its result, stored as attributes of the named
vector which is its primary output. In particular, the `exceedances`
attribute gives upper bounds, not necessarily sharp, for the amount by
which the sum of distances between matched units in the result of
`fullmatch` exceeds the least possible sum of distances between matched
units in a feasible solution to the matching problem given to
`fullmatch`. (Such a bound is also printed by `print.optmatch` and by
`summary.optmatch`.)

## References

Hansen, B.B. and Klopfer, S.O. (2006), ‘Optimal full matching and
related designs via network flows’, *Journal of Computational and
Graphical Statistics*, **15**, 609–627.

## See also

[`matched`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md),
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)

## Examples

``` r
data(nuclearplants)

### Pair matching on a Mahalanobis distance
( pm1 <- pairmatch(pr ~ t1 + t2, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.1  1.1 1.10  1.2 <NA> <NA>  1.3  1.3  1.8 <NA> <NA> <NA> <NA>  1.9 <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7 <NA>  1.6 <NA> <NA> <NA>  1.8  1.9 1.10 
summary(pm1)
#> Structure of matched sets:
#> 1:1 0:1 
#>  10  12 
#> Effective Sample Size:  10 
#> (equivalent number of matched pairs).
#> 

### Pair matching within a propensity score caliper.
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
### For more complicated models, create a distance matrix and pass it to fullmatch.
mhd <- match_on(pr ~ t1 + t2, data = nuclearplants) + caliper(match_on(ppty), 2)
( pm2 <- pairmatch(mhd, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA> 1.10  1.1 <NA>  1.2 <NA> <NA>  1.3  1.3  1.2 <NA> <NA> <NA> <NA>  1.4  1.8 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.9  1.5 <NA>  1.6  1.5  1.7 <NA>  1.6  1.1 <NA> <NA>  1.8  1.9 1.10 
summary(pm2)
#> Structure of matched sets:
#> 1:1 0:1 
#>  10  12 
#> Effective Sample Size:  10 
#> (equivalent number of matched pairs).
#> 

### Propensity balance assessment. Requires RItools package.
if(require(RItools)) summary(pm2, ppty)
#> Structure of matched sets:
#> 1:1 0:1 
#>  10  12 
#> Effective Sample Size:  10 
#> (equivalent number of matched pairs).
#> 
#> Balance test overall result:
#>   chisquare df p.value
#>        8.54  9   0.481

### 1:2 matched triples
( tm <- pairmatch(pr ~ t1 + t2, controls = 2, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#>  1.3  1.1  1.1 1.10  1.2  1.3 1.10  1.5  1.3  1.8  1.8 <NA>  1.4  1.7  1.9  1.1 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7  1.6  1.6  1.9 <NA>  1.2  1.8  1.9 1.10 
summary(tm)
#> Structure of matched sets:
#> 1:2 0:1 
#>  10   2 
#> Effective Sample Size:  13.3 
#> (equivalent number of matched pairs).
#> 

### Creating a data frame with the matched sets attached.
### match_on(), caliper() and the like cooperate with pairmatch()
### to make sure observations are in the proper order:
all.equal(names(tm), row.names(nuclearplants))
#> [1] TRUE
### So our data frame including the matched sets is just
cbind(nuclearplants, matches=tm)
#>     cost  date t1 t2  cap pr ne ct bw cum.n pt matches
#> H 460.05 68.58 14 46  687  0  1  0  0    14  0     1.3
#> I 452.99 67.33 10 73 1065  0  0  1  0     1  0     1.1
#> A 443.22 67.33 10 85 1065  1  0  1  0     1  0     1.1
#> J 652.32 68.00 11 67 1065  0  1  1  0    12  0    1.10
#> B 642.23 68.00 11 78 1065  1  1  1  0    12  0     1.2
#> K 345.39 67.92 13 51  514  0  1  1  0     3  0     1.3
#> L 272.37 68.17 12 50  822  0  0  0  0     5  0    1.10
#> M 317.21 68.42 14 59  457  0  0  0  0     1  0     1.5
#> C 457.12 68.42 15 55  822  1  0  0  0     5  0     1.3
#> N 690.19 68.33 12 71  792  0  1  1  1     2  0     1.8
#> O 350.63 68.58 12 64  560  0  0  0  0     3  0     1.8
#> P 402.59 68.75 13 47  790  0  1  0  0     6  0    <NA>
#> Q 412.18 68.42 15 62  530  0  0  1  0     2  0     1.4
#> R 495.58 68.92 17 52 1050  0  0  0  0     7  0     1.7
#> S 394.36 68.92 13 65  850  0  0  0  1    16  0     1.9
#> T 423.32 68.42 11 67  778  0  0  0  0     3  0     1.1
#> U 712.27 69.50 18 60  845  0  1  0  0    17  0     1.7
#> D 289.66 68.42 15 76  530  1  0  1  0     2  0     1.4
#> V 881.24 69.17 15 67 1090  0  0  0  0     1  0     1.4
#> E 490.88 68.92 16 59 1050  1  0  0  0     8  0     1.5
#> W 567.79 68.75 11 70  913  0  0  1  1    15  0     1.2
#> F 665.99 70.92 22 57  828  1  1  0  0    20  0     1.6
#> X 621.45 69.67 16 59  786  0  0  1  0    18  0     1.5
#> G 608.80 70.08 19 58  821  1  0  0  0     3  0     1.7
#> Y 473.64 70.42 19 44  538  0  0  1  0    19  0     1.6
#> Z 697.14 71.08 20 57 1130  0  0  1  0    21  0     1.6
#> d 207.51 67.25 13 63  745  0  0  0  0     8  1     1.9
#> e 288.48 67.17  9 48  821  0  0  1  0     7  1    <NA>
#> f 284.88 67.83 12 63  886  0  0  0  1    11  1     1.2
#> a 280.36 67.83 12 71  886  1  0  0  1    11  1     1.8
#> b 217.38 67.25 13 72  745  1  0  0  0     8  1     1.9
#> c 270.71 67.83  7 80  886  1  0  0  1    11  1    1.10

### In contrast, if your matching distance is an ordinary matrix
### (as earlier versions of optmatch required), you'll
### have to align it by observation name with your data set.
cbind(nuclearplants, matches = tm[row.names(nuclearplants)])
#>     cost  date t1 t2  cap pr ne ct bw cum.n pt matches
#> H 460.05 68.58 14 46  687  0  1  0  0    14  0     1.3
#> I 452.99 67.33 10 73 1065  0  0  1  0     1  0     1.1
#> A 443.22 67.33 10 85 1065  1  0  1  0     1  0     1.1
#> J 652.32 68.00 11 67 1065  0  1  1  0    12  0    1.10
#> B 642.23 68.00 11 78 1065  1  1  1  0    12  0     1.2
#> K 345.39 67.92 13 51  514  0  1  1  0     3  0     1.3
#> L 272.37 68.17 12 50  822  0  0  0  0     5  0    1.10
#> M 317.21 68.42 14 59  457  0  0  0  0     1  0     1.5
#> C 457.12 68.42 15 55  822  1  0  0  0     5  0     1.3
#> N 690.19 68.33 12 71  792  0  1  1  1     2  0     1.8
#> O 350.63 68.58 12 64  560  0  0  0  0     3  0     1.8
#> P 402.59 68.75 13 47  790  0  1  0  0     6  0    <NA>
#> Q 412.18 68.42 15 62  530  0  0  1  0     2  0     1.4
#> R 495.58 68.92 17 52 1050  0  0  0  0     7  0     1.7
#> S 394.36 68.92 13 65  850  0  0  0  1    16  0     1.9
#> T 423.32 68.42 11 67  778  0  0  0  0     3  0     1.1
#> U 712.27 69.50 18 60  845  0  1  0  0    17  0     1.7
#> D 289.66 68.42 15 76  530  1  0  1  0     2  0     1.4
#> V 881.24 69.17 15 67 1090  0  0  0  0     1  0     1.4
#> E 490.88 68.92 16 59 1050  1  0  0  0     8  0     1.5
#> W 567.79 68.75 11 70  913  0  0  1  1    15  0     1.2
#> F 665.99 70.92 22 57  828  1  1  0  0    20  0     1.6
#> X 621.45 69.67 16 59  786  0  0  1  0    18  0     1.5
#> G 608.80 70.08 19 58  821  1  0  0  0     3  0     1.7
#> Y 473.64 70.42 19 44  538  0  0  1  0    19  0     1.6
#> Z 697.14 71.08 20 57 1130  0  0  1  0    21  0     1.6
#> d 207.51 67.25 13 63  745  0  0  0  0     8  1     1.9
#> e 288.48 67.17  9 48  821  0  0  1  0     7  1    <NA>
#> f 284.88 67.83 12 63  886  0  0  0  1    11  1     1.2
#> a 280.36 67.83 12 71  886  1  0  0  1    11  1     1.8
#> b 217.38 67.25 13 72  745  1  0  0  0     8  1     1.9
#> c 270.71 67.83  7 80  886  1  0  0  1    11  1    1.10


### Match in subgroups only. There are a few ways to specify this.
m1 <- pairmatch(pr ~ t1 + t2, data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m2 <- pairmatch(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
### Matching on propensity scores within matching in subgroups only:
m3 <- pairmatch(glm(pr ~ t1 + t2, data=nuclearplants, family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- pairmatch(glm(pr ~ t1 + t2 + pt, data=nuclearplants,
                    family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- pairmatch(glm(pr ~ t1 + t2 + strata(pt), data=nuclearplants,
                    family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
```
