# Compute value of primal problem given flows and arc costs

Compute value of primal problem given flows and arc costs

## Usage

``` r
evaluate_primal(distances, solution)
```

## Arguments

- distances:

  An InfinitySparseMatrix giving distances

- solution:

  A MCFSolutions object

## Value

The value of the primal problem, i.e. sum of products of `distances`
with flow along arcs in `solution`

## Author

Hansen
