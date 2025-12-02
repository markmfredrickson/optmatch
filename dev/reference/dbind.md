# Diagonally bind together subgroup-specific distances

This function generates a single block-diagonal distance matrix given
several distance matrices defined on subgroups.

## Usage

``` r
dbind(..., force_unique_names = FALSE)
```

## Arguments

- ...:

  Any number of distance objects which can be converted to
  `InfinitySparseMatrix`, such as class `matrix`, `DenseMatrix`,
  `InfinitySparseMatrix`, or `BlockedInfinitySparseMatrix`, or `list`s
  containing distance objects.

- force_unique_names:

  Default `FALSE`. When row or column names are not unique among all
  distances, if `FALSE`, throw a warning and rename all rows and columns
  to ensure unique names. If `TRUE`, error on non-unique names.

## Value

A `BlockedInfinitySparseMatrix` containing a block-diagonal distance
matrix. If only a single distance is passed to `dbind` and it is not
already a `BlockedInfinitySparseMatrix`, the result will be an
`InfinitySparseMatrix` instead.

## Details

When you've generated several distances matrices on subgroups in your
analysis, you may wish to combine them into a single block-diagonal
distance matrix. The `dbind` function facilitates this.

Any `BlockedInfinitySparseMatrix` include in `...` will be broken into
individual `InfinitySparseMatrix` before being joined back together. For
example, if `b` is a `BlockedInfinitySparseMatrix` with 2 subgroups and
`m` is a distance without subgroups, then `dbind(b, m)` will be a
`BlockedInfinitySparseMatrix` with 3 subgroups.

If there are any shared names (either row or column) among all distances
passed in, by default all matrices will be renamed to ensure unique
names by appending "X." to each distance, where "X" is ascending lower
case letters ("a.", "b.", etc). Setting the `force_unique_names`
argument to `TRUE` errors on this instead.

If the matrices need to be renamed and there are more than 26 separate
matrices, after the first 26 single "X." prefixs, they will continue as
"YX.", e.g "aa.", "ab.", "ac.". If more than 676 separate matrices, the
prefix wil continue to "ZYX.", e.g. "aaa.", "aab.", "aac.". This scheme
supports up to 18,278 unique matrices.

Note that you do **not** have to combine subgroup distances into a
single blocked distance using this function to ultimately obtain a
single matching set. Instead, take a look at the vignette
[`vignette("matching-within-subgroups", package = "optmatch")`](https://markmfredrickson.github.io/optmatch/dev/articles/matching-within-subgroups.md)
for details on combining multiple matches.

## Examples

``` r
data(nuclearplants)
m1 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 0),
               caliper = 1)
m2 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 1),
               caliper = 1.3)
blocked <- dbind(m1, m2)

dists <- list(m1, m2)

blocked2 <- dbind(dists)
identical(blocked, blocked2)
#> [1] TRUE
```
