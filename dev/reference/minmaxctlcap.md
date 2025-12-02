# Set thinning and thickening caps for full matching

Functions to find the largest value of min.controls, or the smallest
value of max.controls, for which a full matching problem is feasible.
These are determined by constraints embedded in the matching problem's
distance matrix.

## Usage

``` r
maxControlsCap(distance, min.controls = NULL, solver = "")

minControlsCap(distance, max.controls = NULL, solver = "")
```

## Arguments

- distance:

  Either a matrix of non-negative, numeric discrepancies, or a list of
  such matrices. (See
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  for details.)

- min.controls:

  Optionally, set limits on the minimum number of controls per matched
  set. (Only makes sense for `maxControlsCap`.)

- solver:

  Choose which solver to use. See
  [`help(fullmatch)`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  for details.

- max.controls:

  Optionally, set limits on the maximum number of controls per matched
  set. (Only makes sense for `minControlsCap`.)

## Value

For `minControlsCap`, `strictest.feasible.min.controls` and
`given.max.controls`. For `maxControlsCap`, `given.min.controls` and
`strictest.feasible.max.controls`.

- strictest.feasible.min.controls:

  The largest values of the
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  argument `min.controls` that yield a full match;

- given.max.controls:

  The `max.controls` argument given to `minControlsCap` or, if none was
  given, a vector of `Inf`s.

- given.min.controls:

  The `min.controls` argument given to `maxControlsCap` or, if none was
  given, a vector of `0`s;

- strictest.feasible.max.controls:

  The smallest values of the
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  argument `max.controls` that yield a full match.

## Details

The function works by repeated application of full matching, so on large
problems it can be time-consuming.

## Note

Essentially this is just a line search. I've done several things to
speed it up, but not everything that might be done. At present, not very
thoroughly tested either: you might check the final results to make sure
that
[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
works with the values of `min.controls` (or `max.controls`) suggested by
these functions, and that it ceases to work if you increase (decrease)
those values. Comments appreciated.

## References

Hansen, B.B. and S. Olsen Klopfer (2006), ‘Optimal full matching and
related designs via network flows’, *Journal of Computational and
Graphical Statistics* **15**, 609–627.

## See also

[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)

## Author

Ben B. Hansen
