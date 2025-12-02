# (Internal) Report successful subproblems.

[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
can break up a large matching problem into smaller subproblems (for
example, using strata defined by
[`exactMatch`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)).
This function lists the subproblems in a match and list whether at least
on treated unit was matched in subproblem. Subproblems that have no
matched treated units are said to have "failed."

## Usage

``` r
subproblemSuccess(x)
```

## Arguments

- x:

  The result of
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

## Value

A named logical vector indicating either success or failure for each
subproblem.
