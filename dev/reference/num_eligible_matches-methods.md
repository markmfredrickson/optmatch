# Returns the number of eligible matches for the distance.

This will return a list of the number of finite entries in a distance
matrix. If the distance has no subgroups, it will be a list of length 1.
If the distance has subgroups (i.e. `x` is an
`BlockedInfinitySparseMatrix`, it will be a named list.)

## Usage

``` r
num_eligible_matches(x)

# S3 method for class 'optmatch.dlist'
num_eligible_matches(x)

# S3 method for class 'matrix'
num_eligible_matches(x)

# S3 method for class 'InfinitySparseMatrix'
num_eligible_matches(x)

# S3 method for class 'BlockedInfinitySparseMatrix'
num_eligible_matches(x)
```

## Arguments

- x:

  Any distance object.

## Value

A list counting the number of eligible matches in the distance.
