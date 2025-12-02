# Optimal full matching

Given two groups, such as a treatment and a control group, and a method
of creating a treatment-by-control discrepancy matrix indicating
desirability and permissibility of potential matches (or optionally an
already created such discrepancy matrix), create optimal full matches of
members of the groups. Optionally, incorporate restrictions on matched
sets' ratios of treatment to control units.

## Usage

``` r
fullmatch(
  x,
  min.controls = 0,
  max.controls = Inf,
  omit.fraction = NULL,
  mean.controls = NULL,
  tol = 0.001,
  data = NULL,
  solver = "",
  ...
)

full(
  x,
  min.controls = 0,
  max.controls = Inf,
  omit.fraction = NULL,
  mean.controls = NULL,
  tol = 0.001,
  data = NULL,
  solver = "",
  ...
)
```

## Arguments

- x:

  Any valid input to `match_on`. `fullmatch` will use `x` and any
  optional arguments to generate a distance before performing the
  matching.

  If `x` is a numeric vector, there must also be passed a vector `z`
  indicating grouping. Both vectors must be named.

  Alternatively, a precomputed distance may be entered. A matrix of
  non-negative discrepancies, each indicating the permissibility and
  desirability of matching the unit corresponding to its row (a
  'treatment') to the unit corresponding to its column (a 'control');
  or, better, a distance specification as produced by
  [`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md).

- min.controls:

  The minimum ratio of controls to treatments that is to be permitted
  within a matched set: should be non-negative and finite. If
  `min.controls` is not a whole number, the reciprocal of a whole
  number, or zero, then it is rounded *down* to the nearest whole number
  or reciprocal of a whole number.

  When matching within subclasses (such as those created by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)),
  `min.controls` may be a named numeric vector separately specifying the
  minimum permissible ratio of controls to treatments for each subclass.
  The names of this vector should include names of all subproblems
  `distance`.

- max.controls:

  The maximum ratio of controls to treatments that is to be permitted
  within a matched set: should be positive and numeric. If
  `max.controls` is not a whole number, the reciprocal of a whole
  number, or `Inf`, then it is rounded *up* to the nearest whole number
  or reciprocal of a whole number.

  When matching within subclasses (such as those created by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)),
  `max.controls` may be a named numeric vector separately specifying the
  maximum permissible ratio of controls to treatments in each subclass.

- omit.fraction:

  Optionally, specify what fraction of controls or treated subjects are
  to be rejected. If `omit.fraction` is a positive fraction less than
  one, then `fullmatch` leaves up to that fraction of the control
  reservoir unmatched. If `omit.fraction` is a negative number greater
  than -1, then `fullmatch` leaves up to \|`omit.fraction`\| of the
  treated group unmatched. Positive values are only accepted if
  `max.controls` \>= 1; negative values, only if `min.controls` \<= 1.
  If neither `omit.fraction` or `mean.controls` are specified, then only
  those treated and control subjects without permissible matches among
  the control and treated subjects, respectively, are omitted.

  When matching within subclasses (such as those created by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)),
  `omit.fraction` specifies the fraction of controls to be rejected in
  each subproblem, a parameter that can be made to differ by subclass by
  setting `omit.fraction` equal to a named numeric vector of fractions.

  At most one of `mean.controls` and `omit.fraction` can be non-`NULL`.

- mean.controls:

  Optionally, specify the average number of controls per treatment to be
  matched. Must be no less than than `min.controls` and no greater than
  the either `max.controls` or the ratio of total number of controls
  versus total number of treated. Some controls will likely not be
  matched to ensure meeting this value. If neither `omit.fraction` or
  `mean.controls` are specified, then only those treated and control
  subjects without permissible matches among the control and treated
  subjects, respectively, are omitted.

  When matching within subclasses (such as those created by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)),
  `mean.controls` specifies the average number of controls per treatment
  per subproblem, a parameter that can be made to differ by subclass by
  setting `mean.controls` equal to a named numeric vector.

  At most one of `mean.controls` and `omit.fraction` can be non-`NULL`.

- tol:

  Because of internal rounding, `fullmatch` may solve a slightly
  different matching problem than the one specified, in which the match
  generated by `fullmatch` may not coincide with an optimal solution of
  the specified problem. `tol` times the number of subjects to be
  matched specifies the extent to which `fullmatch`'s output is
  permitted to differ from an optimal solution to the original problem,
  as measured by the sum of discrepancies for all treatments and
  controls placed into the same matched sets.

- data:

  Optional `data.frame` or `vector` to use to get order of the final
  matching factor. If a `data.frame`, the `rownames` are used. If a
  vector, the `names` are first tried, otherwise the contents is
  considered to be a character vector of names. Useful to pass if you
  want to combine a match (using, e.g., `cbind`) with the data that were
  used to generate it (for example, in a propensity score matching).

- solver:

  Choose which solver to use. Currently implemented are RELAX-IV and
  LEMON. Default of `""`, a blank string, will use RELAX-IV if the
  **rrelaxiv** package is installed, otherwise will use LEMON.

  To explicitly use RELAX-IV, pass string "RELAX-IV".

  To use LEMON, pass string "LEMON". Optionally, to specify which
  algorithm LEMON will use, pass the function
  [LEMON](https://markmfredrickson.github.io/optmatch/dev/reference/LEMON.md)
  with argument for the algorithm name, "CycleCancelling",
  "CapacityScaling", "CostScaling", and "NetworkSimplex". See this site
  for details on their differences:
  <https://lemon.cs.elte.hu/pub/doc/latest/a00606.html>. CycleCancelling
  is the default.

  The CycleCancelling algorithm seems to produce results most closely
  resembling those of optmatch versions prior to 1.0. We have observed
  the other LEMON algorithms to produce different results when the
  `mean.controls` is unspecified, or specified in such a way as to
  produce an infeasible matching problem. When using a LEMON algorithm
  other than CycleCancelling, we recommend setting the
  `fullmatch_try_recovery` option to `FALSE`.

- ...:

  Additional arguments, passed to `match_on` (e.g. `within`) or to
  specific methods.

## Value

A
[`optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
object (`factor`) indicating matched groups.

## Details

If passing an already created discrepancy matrix, finite entries
indicate permissible matches, with smaller discrepancies indicating more
desirable matches. The matrix must have row and column names.

If it is desirable to create the discrepancies matrix beforehand (for
example, if planning on running several different matching schemes),
consider using
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
to generate the distances. This generic function has several useful
methods for handling propensity score models, computing Mahalanobis
distances (and other arbitrary distances), and using user supplied
functions. These distances can also be combined with those generated by
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
and
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
to create very nuanced matching specifications.

The value of `tol` can have a substantial effect on computation time;
with smaller values, computation takes longer. Not every tolerance can
be met, and how small a tolerance is too small varies with the machine
and with the details of the problem. If `fullmatch` can't guarantee that
the tolerance is as small as the given value of argument `tol`, then
matching proceeds but a warning is issued.

By default, `fullmatch` will attempt, if the given constraints are
infeasible, to find a feasible problem using the same constraints. This
will almost surely involve using a more restrictive `omit.fraction` or
`mean.controls`. (This will never automatically omit treatment units.)
Note that this does not guarantee that the returned match has the least
possible number of omitted subjects, it only gives a match that is
feasible within the given constraints. It may often be possible to
loosen the `omit.fraction` or `mean.controls` constraint and still find
a feasible match. The auto recovery is controlled by
`options("fullmatch_try_recovery")`.

In full matching problems permitting many-one matches (`min.controls`
less than 1), the number of controls contributing to matches can exceed
what was requested by setting a value of `mean.controls` or
`omit.fraction`. I.e., in this setting `mean.controls` sets the minimum
ratio of number of controls to number of treatments placed into matched
sets.

If the program detects that (what it thinks is) a large problem, a
warning is issued. Unless you have an older computer, there's a good
chance that you can handle larger problems (at the cost of increased
computation time). To check the large problem threshold, use
[`getMaxProblemSize`](https://markmfredrickson.github.io/optmatch/dev/reference/getMaxProblemSize.md);
to re-set it, use
[`setMaxProblemSize`](https://markmfredrickson.github.io/optmatch/dev/reference/setMaxProblemSize.md).

## References

Hansen, B.B. and Klopfer, S.O. (2006), ‘ Optimal full matching and
related designs via network flows’, *Journal of Computational and
Graphical Statistics*, **15**, 609–627.

Hansen, B.B. (2004), ‘Full Matching in an Observational Study of
Coaching for the SAT’, *Journal of the American Statistical
Association*, **99**, 609–618.

Rosenbaum, P. (1991), ‘A Characterization of Optimal Designs for
Observational Studies’, *Journal of the Royal Statistical Society,
Series B*, **53**, 597–610.

## Examples

``` r
data(nuclearplants)
### Full matching on a Mahalanobis distance.
( fm1 <- fullmatch(pr ~ t1 + t2, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#>  1.3  1.1  1.1  1.8  1.2  1.3  1.3  1.3  1.3  1.8  1.8  1.3  1.5  1.3  1.9 1.10 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7  1.3  1.6  1.3  1.3  1.3  1.8  1.9 1.10 
summary(fm1)
#> Structure of matched sets:
#>  1:1  1:2  1:3 1:5+ 
#>    7    1    1    1 
#> Effective Sample Size:  11.7 
#> (equivalent number of matched pairs).
#> 

### Full matching with restrictions.
( fm2 <- fullmatch(pr ~ t1 + t2, min.controls = .5, max.controls = 4, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#>  1.3  1.1  1.1  1.8  1.2  1.3  1.3  1.5  1.3  1.8  1.8  1.3  1.5  1.5  1.9  1.8 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7  1.7  1.6  1.9 1.10  1.9  1.8  1.9 1.10 
summary(fm2)
#> Structure of matched sets:
#> 1:1 1:2 1:3 1:4 
#>   5   1   1   3 
#> Effective Sample Size:  12.6 
#> (equivalent number of matched pairs).
#> 

### Full matching to half of available controls.
( fm3 <- fullmatch(pr ~ t1 + t2, omit.fraction = .5, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.1  1.1 1.10  1.2 <NA> <NA>  1.3  1.3  1.8 <NA> <NA>  1.5 <NA>  1.9 <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7 <NA>  1.6 <NA> <NA> <NA>  1.8  1.9 1.10 
summary(fm3)
#> Structure of matched sets:
#> 1:1 1:2 0:1 
#>   9   1  11 
#> Effective Sample Size:  10.3 
#> (equivalent number of matched pairs).
#> 

### Full matching attempts recovery when the initial restrictions are infeasible.
### Limiting max.controls = 1 allows use of only 10 of 22 controls.
( fm4 <- fullmatch(pr ~ t1 + t2, max.controls = 1, data=nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.1  1.1 1.10  1.2 <NA> <NA>  1.3  1.3  1.8 <NA> <NA> <NA> <NA>  1.9 <NA> 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7 <NA>  1.6 <NA> <NA> <NA>  1.8  1.9 1.10 
summary(fm4)
#> Structure of matched sets:
#> 1:1 0:1 
#>  10  12 
#> Effective Sample Size:  10 
#> (equivalent number of matched pairs).
#> 
### To recover restrictions
optmatch_restrictions(fm4)
#> $min.controls
#> [1] 0
#> 
#> $max.controls
#> [1] 1
#> 
#> $omit.fraction
#> [1] 0.5454545
#> 

### Full matching within a propensity score caliper.
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
### Note that units without counterparts within the caliper are automatically dropped.
### For more complicated models, create a distance matrix and pass it to fullmatch.
mhd <- match_on(pr ~ t1 + t2, data = nuclearplants) + caliper(match_on(ppty), width = 1)
( fm5 <- fullmatch(mhd, data = nuclearplants) )
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#> <NA>  1.9 <NA>  1.2  1.1 <NA> <NA>  1.2  1.2  1.2  1.2 <NA>  1.2  1.2  1.4  1.1 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.5  1.3  1.3  1.4  1.2  1.5  1.2  1.5 <NA>  1.4  1.7 <NA>  1.2  1.7  1.3  1.9 
summary(fm5)
#> Structure of matched sets:
#>  1:0  2:1  1:1  1:2 1:5+  0:1 
#>    1    2    3    1    1    6 
#> Effective Sample Size:  8.8 
#> (equivalent number of matched pairs).
#> 

### Propensity balance assessment. Requires RItools package.
if (require(RItools)) summary(fm5,ppty)
#> Loading required package: RItools
#> Loading required package: ggplot2
#> Loading required package: survival
#> 
#> Attaching package: ‘survival’
#> The following object is masked from ‘package:optmatch’:
#> 
#>     strata
#> Structure of matched sets:
#>  1:0  2:1  1:1  1:2 1:5+  0:1 
#>    1    2    3    1    1    6 
#> Effective Sample Size:  8.8 
#> (equivalent number of matched pairs).
#> 
#> Balance test overall result:
#>   chisquare df p.value
#>        7.93  9   0.542

### The order of the names in the match factor is the same
### as the nuclearplants data.frame since we used the data argument
### when calling fullmatch. The order would be unspecified otherwise.
cbind(nuclearplants, matches = fm5)
#>     cost  date t1 t2  cap pr ne ct bw cum.n pt matches
#> H 460.05 68.58 14 46  687  0  1  0  0    14  0    <NA>
#> I 452.99 67.33 10 73 1065  0  0  1  0     1  0     1.9
#> A 443.22 67.33 10 85 1065  1  0  1  0     1  0    <NA>
#> J 652.32 68.00 11 67 1065  0  1  1  0    12  0     1.2
#> B 642.23 68.00 11 78 1065  1  1  1  0    12  0     1.1
#> K 345.39 67.92 13 51  514  0  1  1  0     3  0    <NA>
#> L 272.37 68.17 12 50  822  0  0  0  0     5  0    <NA>
#> M 317.21 68.42 14 59  457  0  0  0  0     1  0     1.2
#> C 457.12 68.42 15 55  822  1  0  0  0     5  0     1.2
#> N 690.19 68.33 12 71  792  0  1  1  1     2  0     1.2
#> O 350.63 68.58 12 64  560  0  0  0  0     3  0     1.2
#> P 402.59 68.75 13 47  790  0  1  0  0     6  0    <NA>
#> Q 412.18 68.42 15 62  530  0  0  1  0     2  0     1.2
#> R 495.58 68.92 17 52 1050  0  0  0  0     7  0     1.2
#> S 394.36 68.92 13 65  850  0  0  0  1    16  0     1.4
#> T 423.32 68.42 11 67  778  0  0  0  0     3  0     1.1
#> U 712.27 69.50 18 60  845  0  1  0  0    17  0     1.5
#> D 289.66 68.42 15 76  530  1  0  1  0     2  0     1.3
#> V 881.24 69.17 15 67 1090  0  0  0  0     1  0     1.3
#> E 490.88 68.92 16 59 1050  1  0  0  0     8  0     1.4
#> W 567.79 68.75 11 70  913  0  0  1  1    15  0     1.2
#> F 665.99 70.92 22 57  828  1  1  0  0    20  0     1.5
#> X 621.45 69.67 16 59  786  0  0  1  0    18  0     1.2
#> G 608.80 70.08 19 58  821  1  0  0  0     3  0     1.5
#> Y 473.64 70.42 19 44  538  0  0  1  0    19  0    <NA>
#> Z 697.14 71.08 20 57 1130  0  0  1  0    21  0     1.4
#> d 207.51 67.25 13 63  745  0  0  0  0     8  1     1.7
#> e 288.48 67.17  9 48  821  0  0  1  0     7  1    <NA>
#> f 284.88 67.83 12 63  886  0  0  0  1    11  1     1.2
#> a 280.36 67.83 12 71  886  1  0  0  1    11  1     1.7
#> b 217.38 67.25 13 72  745  1  0  0  0     8  1     1.3
#> c 270.71 67.83  7 80  886  1  0  0  1    11  1     1.9

### Match in subgroups only. There are a few ways to specify this.
m1 <- fullmatch(pr ~ t1 + t2, data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m2 <- fullmatch(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
### Matching on propensity scores within matching in subgroups only:
m3 <- fullmatch(glm(pr ~ t1 + t2, data=nuclearplants, family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- fullmatch(glm(pr ~ t1 + t2 + pt, data=nuclearplants,
                    family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- fullmatch(glm(pr ~ t1 + t2 + strata(pt), data=nuclearplants,
                    family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
```
