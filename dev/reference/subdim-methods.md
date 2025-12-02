# Returns the dimension of each valid subproblem

Returns a list containing the dimensions of all valid subproblems.

## Usage

``` r
subdim(x)

# S3 method for class 'InfinitySparseMatrix'
subdim(x)

# S3 method for class 'matrix'
subdim(x)

# S3 method for class 'BlockedInfinitySparseMatrix'
subdim(x)

# S3 method for class 'optmatch.dlist'
subdim(x)
```

## Arguments

- x:

  A distance specification to get the sub-dimensions of.

## Value

A data frame listing the dimensions of each valid subproblem. Any
subproblems with 0 controls or 0 treatments will be ignored. The names
of the entries in the list will be the names of the subproblems, if they
exist. There will be two rows, named "treatments" and "controls".

## Examples

``` r
em <- exactMatch(pr ~ pt, data=nuclearplants)
m1 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants)
stratumStructure(m1)
#>  1:1  1:2  1:3 1:10 
#>    7    1    1    1 
(subdims_em <- subdim(em))
#>             0 1
#> treatments  7 3
#> controls   19 3
m2 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants,
                mean.controls=pmin(1.5, subdims_em["controls",] / subdims_em["treatments",])
                )
stratumStructure(m2)
#> 1:1 1:2 1:3 0:1 
#>   7   2   1   8 
                
```
