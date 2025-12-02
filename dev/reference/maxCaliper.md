# Find the maximum caliper width that will create a feasible problem.

Larger calipers permit more possible matches between treated and control
groups, which can be better for creating matches with larger effective
sample sizes. The downside is that wide calipers may make the matching
problem too big for processor or memory constraints. `maxCaliper`
attempts to find a caliper value, for a given vector of scores and a
treatment indicator, that will be possible given the maximum problem
size constraints imposed by
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
and
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

## Usage

``` r
maxCaliper(scores, z, widths, structure = NULL, exact = TRUE)
```

## Arguments

- scores:

  A numeric vector of scores providing 1-D position of units

- z:

  Treatment indicator vector

- widths:

  A vector of caliper widths to try, will be sorted largest to smallest.

- structure:

  Optional factor variable that groups the scores, as would be used by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md).
  Including structure allows for wider calipers.

- exact:

  A logical indicating if the exact problem size should be computed
  (`exact = TRUE`) or if a more computationally efficient upper bound
  should be used instead (`exact = FALSE`). The upper bound may lead to
  narrower calipers, even if wider calipers would have sufficed using
  the exact method.

## Value

numeric The value of the largest caliper that creates a feasible
problem. If no such caliper exists in `widths`, an error will be
generated.
