# (Internal) Create `optmatch` objects, the result of matching.

This internal function is used to create the final output of the
matching functions
([`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
and
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)).
The `optmatch` object descends from a `factor`, but contains additional
information relating to the quality of the match.

## Usage

``` r
makeOptmatch(distance, solutions, call, data = NULL)
```

## Arguments

- distance:

  A `DistanceSpecificaton` object used to create the match.

- solutions:

  A list of the results of the matching, one `list(cells,maxErr)` object
  per subproblem.

- call:

  The call to `fullmatch` or `pairmatch` to be displayed later.

- data:

  An object from which `names` or `row.names` will provide the order of
  the items in the match. If no names are attached to this object, the
  contents will be used as names.

## Value

`optmatch` object

## See also

[`summary.optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
