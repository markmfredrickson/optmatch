# Set the maximum problem size

Helper function to ease setting the largest problem size to be accepted
by `pairmatch` or `fullmatch`.

## Usage

``` r
setMaxProblemSize(size = Inf)
```

## Arguments

- size:

  Positive integer, or `Inf`

## Details

The function sets the optmatch_max_problem_size global option. The
option ships with the option pre-set to a value that's relatively small,
smaller than what most modern computers can handle. Invoking this
function with no argument re-sets the optmatch_max_problem_size option
to `Inf`, effectively disabling checks on problem size. Unless you're
working with an older computer, it probably makes sense for most users
to do this, at least until they determine what problem sizes are too
large for their machines. (You'll know that when your R crashes, or
simply takes too long for your taste.)

To determine the size of a problem without subproblems, i.e. exact
matching categories, use
[`match_on`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
to set up and store the problem distance, then apply `length` to the
result. If there were exact matching constraints imposed during the
creation of the distance, then you'll want to look at the largest size
of a subproblem. Apply
[`findSubproblems`](https://markmfredrickson.github.io/optmatch/dev/reference/findSubproblems.md)
to your distance, creating a list, say `dlist`, of your distances; then
do `sapply(dlist, length)` to determine the sizes of the subproblems.

## See also

[`getMaxProblemSize`](https://markmfredrickson.github.io/optmatch/dev/reference/getMaxProblemSize.md)

## Author

Ben B. Hansen
