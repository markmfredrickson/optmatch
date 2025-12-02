# Checks if two distances are equivalent. `x` and `y` can be distances (`InfinitySparseMatrix`, `BlockedInfinitySparseMatrix`, or `DenseMatrix`), or they can be `optmatch` objects.

To save space, `optmatch` objects merely store a hash of the distance
matrix instead of the original object. Any distance objects are hashed
before comparison.

## Usage

``` r
optmatch_same_distance(x, y)
```

## Arguments

- x:

  A distances (`InfinitySparseMatrix`, `BlockedInfinitySparseMatrix`, or
  `DenseMatrix`), or `optmatch` object.

- y:

  A distances (`InfinitySparseMatrix`, `BlockedInfinitySparseMatrix`, or
  `DenseMatrix`), or `optmatch` object.

## Value

Boolean whether the two distance specifications are identical.

## Details

Note that the distance is hashed with its `call` set to `NULL`. (This
avoids issues where, for example, `match_on(Z~X, data=d, caliper=NULL)`
and `match_on(Z~X, data=d)` produce identical matches but have differing
`call`s.)
