# Blocked Infinity Sparse Matrix

Blocked Infinity Sparse Matrices are similar to Infinity Sparse
Matrices, but they also keep track of the groups of units via an
additional slot, `groups`

## Slots

- `groups`:

  factor vector containing groups, with unit names as labels, when
  possible

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
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`InfinitySparseMatrix-class`](https://markmfredrickson.github.io/optmatch/dev/reference/InfinitySparseMatrix-class.md)

## Author

Mark M. Fredrickson
