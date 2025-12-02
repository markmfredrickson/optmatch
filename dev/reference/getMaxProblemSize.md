# What is the maximum allowed problem size?

To prevent users from starting excessively large matching problems, the
maximum problem size is limited by
`options("optmatch_max_problem_size")`. This function a quick helper to
assist fetching this value as a scalar. If the option isn't set, the
function falls back to the default value, hard coded in the `optmatch`
package.

## Usage

``` r
getMaxProblemSize()
```

## Value

logical

## See also

[`options`](https://rdrr.io/r/base/options.html),
[`setMaxProblemSize`](https://markmfredrickson.github.io/optmatch/dev/reference/setMaxProblemSize.md)

## Examples

``` r
optmatch:::getMaxProblemSize() > 1 & optmatch:::getMaxProblemSize() < 1e100
#> [1] TRUE
```
