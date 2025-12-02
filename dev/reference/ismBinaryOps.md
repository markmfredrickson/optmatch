# Element-wise addition

`e1 + e2` returns the element-wise sum of two InfinitySparseMatrix
objects. If either element is inf then the resulting element will be
inf.

`e1 - e2` returns the element-wise subtraction of two
InfinitySparseMatrix objects. If either element is inf then the
resulting element will be inf.

`e1 * e2` returns the element-wise multiplication of two
InfinitySparseMatrix objects. If either element is inf then the
resulting element will be inf.

`e1 / e2` returns the element-wise division of two InfinitySparseMatrix
objects. If either element is inf then the resulting element will be
inf.

## Usage

``` r
# S4 method for class 'InfinitySparseMatrix,InfinitySparseMatrix'
e1 + e2

# S4 method for class 'InfinitySparseMatrix,InfinitySparseMatrix'
e1 - e2

# S4 method for class 'InfinitySparseMatrix,InfinitySparseMatrix'
e1 * e2

# S4 method for class 'InfinitySparseMatrix,InfinitySparseMatrix'
e1/e2
```

## Arguments

- e1:

  an InfinitySparseMatrix object

- e2:

  an InfinitySparseMatrix object

## Value

an InfinitySparseMatrix object representing the element-wise sum of the
two ISM objects
