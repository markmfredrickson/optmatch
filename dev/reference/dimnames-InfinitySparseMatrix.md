# Get and set dimnames for InfinitySparseMatrix objects

InfinitySparseMatrix objects represent sparse matching problems with
treated units as rows of a matrix and controls units as the columns of
the matrix. The names of the units can be retrieved and set using these
methods.

## Usage

``` r
# S4 method for class 'InfinitySparseMatrix'
dimnames(x)

# S4 method for class 'InfinitySparseMatrix,list'
dimnames(x) <- value

# S4 method for class 'InfinitySparseMatrix,NULL'
dimnames(x) <- value
```

## Arguments

- x:

  An InfinitySparseMatrix object.

- value:

  A list with two entries: the treated names and control names,
  respectively.

## Value

A list with treated and control names.
