# Reset implicit node labels of an ArcInfo object

Reset implicit node labels of an ArcInfo object

## Usage

``` r
revise_ArcInfo_nodelabels(x, new, old_positions = 1L:length(new))
```

## Arguments

- x:

  an ArcInfo object

- new:

  character; the new node labels (level sets for factors encoding arc
  start or end nodes)

- old_positions:

  integer; positions for old levels with new levels vector

## Value

ArcInfo object with new levels

## Author

Ben Hansen
