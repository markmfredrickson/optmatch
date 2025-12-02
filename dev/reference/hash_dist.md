# (Internal) Hashing functions for various distance objects \#

(Internal) Hashing functions for various distance objects \#

## Usage

``` r
hash_dist(dist)

# S3 method for class 'DenseMatrix'
hash_dist(dist)

# S3 method for class 'matrix'
hash_dist(dist)

# S3 method for class 'InfinitySparseMatrix'
hash_dist(dist)

# S3 method for class 'BlockedInfinitySparseMatrix'
hash_dist(dist)

# S3 method for class 'optmatch.dlist'
hash_dist(dist)
```

## Arguments

- dist:

  Distance object to hash. Must be one of `InfinitySparseMatrix`,
  `BlockedInfinitySparseMatrix`, `DenseMatrix`, `matrix`, or
  `distmatch.dlist`.

## Value

Hash on the distance object with a null `call`
