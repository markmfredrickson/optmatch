# (Internal) Returns a reasonable upper bound on the arcs remaining after placing a caliper.

(Internal) Returns a reasonable upper bound on the arcs remaining after
placing a caliper.

## Usage

``` r
caliperUpperBound(scores, z, width, structure = NULL)
```

## Arguments

- scores:

  A numeric vector of scores providing 1-D position of units

- z:

  Treatment indicator vector

- width:

  Width of caliper, must be positive.

- structure:

  Optional factor variable that groups the scores, as would be used by
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md).
  Including structure allows for wider calipers.

## Value

numeric Total number of pairwise distances remaining after the caliper
is placed.
