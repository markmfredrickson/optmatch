# Combine multiple distance specifications into a single distance specification.

Creates a new distance specification from the union of two or more
distance specifications. The constituent distances specifications may
have overlapping treated and control units (identified by the `rownames`
and `colnames` respectively).

## Usage

``` r
distUnion(...)
```

## Arguments

- ...:

  The distance specifications (as created with with
  [`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  [`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
  or other distance creation function).

## Value

An InfinitySparseMatrix object with all treated and control units from
the arguments combined. Duplicate entries are resolved in favor of the
earliest argument (e.g., `distUnion(A, B)` will favor entries in `A`
over entries in `B`).

## Details

For combining multiple distance specifications with common controls, but
different treated units, [`rbind`](https://rdrr.io/r/base/cbind.html)
provides a way to combine the different objects. Likewise,
[`cbind`](https://rdrr.io/r/base/cbind.html) provides a way to combine
distance specifications over common treated units, but different control
units.

`distUnion` can combine distance units that have common treated and
control units into a coherent single distance object. If there are
duplicate treated-control entries in multiple input distances, the first
entry will be used.

## See also

[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md),
[`cbind`](https://rdrr.io/r/base/cbind.html),
[`rbind`](https://rdrr.io/r/base/cbind.html)
