# Objects for sparse matching problems.

`InfinitySparseMatrix` is a special class of distance specifications.
Finite entries indicate possible matches, while infinite or NA entries
indicated non-allowed matches. This data type can be more space
efficient for sparse matching problems. Usually, users will create
distance specification using
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
or
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md).
The ordering of units in an `InfinitySparseMatrix` is not guaranteed to
be maintained after subsetting and/or other operations are performed.

## Slots

- `colnames`:

  vector containing names for all control units. This will either be a
  character vector or NULL if units have no names

- `rownames`:

  vector containing names for all treated units. This will either be a
  character vector or NULL if units have no names

- `cols`:

  vector of integers corresponding to control units

- `rows`:

  vector of integers corresponding to treated units

- `dimension`:

  integer vector containing the number of treated and control units, in
  that order

- `call`:

  function call used to create the `InfinitySparseMatrix`

## See also

[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`caliper`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

## Author

Mark M. Fredrickson
