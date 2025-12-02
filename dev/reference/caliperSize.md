# (Internal) Determines how many other units fall within a caliper distance

The matching functions
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
and
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
have a maximum problem size, based on the number of comparisons between
treated and control units. For a completely dense problem, in which
every treated units is compared to every control unit there are
`length(treated) * length(control)` comparisons. A caliper restricts
which comparisons are valid, disallowing matches of treated and control
pairs that are too far apart. A caliper can significantly decrease the
size of a matching problem. The `caliperSize` function reports exactly
who many valid treated-control comparisons remain after applying a
caliper of the given width.

## Usage

``` r
caliperSize(scores, z, width, structure = NULL)
```

## Arguments

- scores:

  A numeric vector of scores providing 1-D position of units

- z:

  Treatment indicator vector

- width:

  Width of caliper, must be positive

- structure:

  Grouping factor to use in computation

## Value

numeric Total number of pairwise distances remaining after the caliper
is placed.
