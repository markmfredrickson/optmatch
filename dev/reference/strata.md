# Identify Stratafication Variables

This is a special function used only in identifying the strata variables
when defining an `exactMatch` during a call to `fullmatch`, `pairmatch`,
or `match_on`. It should not be called externally.

## Usage

``` r
strata(...)
```

## Arguments

- ...:

  any number of variables of the same length.

## Value

the variables with appropriate labels

## Examples

``` r
data(nuclearplants)
fullmatch(pr ~ cost + strata(pt), data = nuclearplants)
#>   H   I   A   J   B   K   L   M   C   N   O   P   Q   R   S   T   U   D   V   E 
#> 0.3 0.3 0.1 0.2 0.2 0.4 0.4 0.4 0.3 0.6 0.4 0.1 0.1 0.5 0.1 0.1 0.6 0.4 0.6 0.5 
#>   W   F   X   G   Y   Z   d   e   f   a   b   c 
#> 0.7 0.6 0.7 0.7 0.3 0.6 1.2 1.1 1.3 1.1 1.2 1.3 
```
