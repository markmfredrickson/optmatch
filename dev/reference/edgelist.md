# Create EdgeList object

If y is a named character vector, then the names should correspond to
whatever in x would otherwise (i.e. if y were NULL) translate to the
levels set of the nodes-representing columns, while the values
themselves give the new levels.

## Usage

``` r
edgelist(x, y = NULL)
```

## Arguments

- x:

  object to convert to edgelist

- y:

  named character vector giving levels for nodes-representing columns,
  or NULL

## Value

EdgeList

## Author

Ben Hansen
