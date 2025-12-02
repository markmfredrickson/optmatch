# Determine distances between matched units

From a match (as produced by `pairmatch` or `fullmatch`) and a distance,
extract the distances of matched units from their matched counterparts.

## Usage

``` r
matched.distances(matchobj, distance, preserve.unit.names = FALSE)
```

## Arguments

- matchobj:

  Value of a call to `pairmatch` or `fullmatch`.

- distance:

  Either a distance matrix or the value of a call to or `match_on`.

- preserve.unit.names:

  Logical. If TRUE, for each matched set `matched.distances` returns the
  submatrix of the distance matrix corresponding to it; if FALSE, a
  vector containing the distances in that submatrix is returned.

## Value

A list of numeric vectors (or matrices) of distances, one for each
matched set. Note that a matched set with 1 treatment and k controls, or
with k treatments and 1 control, has k, not k+1, distances.

## Details

From a match (as produced by `pairmatch` or `fullmatch`) and a distance,
extract the distances of matched units from their matched counterparts.

## Author

Ben B. Hansen

## Examples

``` r
data(plantdist)
plantsfm <- fullmatch(plantdist)
#> Warning: Without 'data' argument the order of the match is not guaranteed
#>     to be the same as your original data.
(plantsfm.d <- matched.distances(plantsfm,plantdist,pres=TRUE))
#> $`1.1`
#> I 
#> 0 
#> 
#> $`1.2`
#> J 
#> 0 
#> 
#> $`1.3`
#> L N P 
#> 4 6 9 
#> 
#> $`1.4`
#> H K M O Q T 
#> 7 8 2 6 0 4 
#> 
#> $`1.5`
#> R S V W 
#> 0 2 8 4 
#> 
#> $`1.6`
#>  U  Z 
#>  5 12 
#> 
#> $`1.7`
#> X Y 
#> 4 8 
#> 
unlist(lapply(plantsfm.d, max))
#> 1.1 1.2 1.3 1.4 1.5 1.6 1.7 
#>   0   0   9   8   8  12   8 
mean(unlist(plantsfm.d))
#> [1] 4.684211
```
