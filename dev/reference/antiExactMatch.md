# Specify a matching problem where units in a common factor cannot be matched.

This function builds a distance specification where treated units are
infinitely far away from control units that share the same level of a
given factor variable. This can be useful for ensuring that matched
groups come from qualitatively different groups.

## Usage

``` r
antiExactMatch(x, z)
```

## Arguments

- x:

  A factor across which matches should be allowed.

- z:

  A logical or binary vector the same length as `x` indicating treatment
  and control for each unit in the study. TRUE or 1 represents a
  treatment unit, FALSE of 0 represents a control unit. NA units are
  excluded.

## Value

A distance specification that encodes the across factor level
constraint.

## Details

The
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
function provides a way of specifying a matching problem where only
units within a factor level may be matched. This function provides the
reverse scenario: a matching problem in which only units across factor
levels are permitted to match. Like
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
the results of this function will most often be used as a `within`
argument to
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
or another distance specification creation function to limit the scope
of the final distance specification (i.e., disallowing any match between
units with the same value on the factor variable `x`).

## See also

[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

## Examples

``` r
data(nuclearplants)

# force entries to be within the same factor:
em <- fullmatch(exactMatch(pr ~ pt, data = nuclearplants), data = nuclearplants)
table(nuclearplants$pt, em)
#>    em
#>     0.1 0.2 0.3 0.4 0.5 0.6 0.7 1.1 1.2 1.3
#>   0   2   2   2   2   2   2  14   0   0   0
#>   1   0   0   0   0   0   0   0   2   2   2

# force treated and control units to have different values of `pt`:
z <- nuclearplants$pr
names(z) <- rownames(nuclearplants)
aem <- fullmatch(antiExactMatch(nuclearplants$pt, z), data = nuclearplants)
table(nuclearplants$pt, aem)
#>    aem
#>     1.1 1.10 1.2 1.3 1.8 1.9
#>   0   4   16   1   2   2   1
#>   1   1    1   1   1   1   1
```
