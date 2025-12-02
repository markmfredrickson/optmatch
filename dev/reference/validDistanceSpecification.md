# (Internal) Validate that objects are valid distance specifications.

The functions
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
and
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
create optimal matches of treated and control units given a matrix (or
similar representation) of distances between treated and control units.
These distance specifications must implement certain generic functions.
This function checks that all necessary methods exist and the object can
be used to specify distances in a way that the matching functions can
use.

## Usage

``` r
validDistanceSpecification(distance, stopOnProblem = TRUE)
```

## Arguments

- distance:

  The object to test.

- stopOnProblem:

  If `TRUE` (default) the function will raise an error for invalid
  objects. Otherwise, it returns a logical.

## Value

logical
