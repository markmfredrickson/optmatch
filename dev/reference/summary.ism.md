# Summarize a distance matrix

Given a distance matrix, return information above it, including
dimension, sparsity information, unmatchable members, summary of finite
distances, and, in the case of `BlockedInfinitySparseMatrix`, block
structure.

## Usage

``` r
# S3 method for class 'InfinitySparseMatrix'
summary(object, ..., distanceSummary = TRUE)

# S3 method for class 'BlockedInfinitySparseMatrix'
summary(
  object,
  ...,
  distanceSummary = TRUE,
  printAllBlocks = FALSE,
  blockStructure = TRUE
)

# S3 method for class 'DenseMatrix'
summary(object, ..., distanceSummary = TRUE)
```

## Arguments

- object:

  A `InfinitySparseMatrix`, `BlockedInfinitySparseMatrix` or
  `DenseMatrix`.

- ...:

  Ignored.

- distanceSummary:

  Default `TRUE`. Should a summary of minimum distance per treatment
  member be calculated? May be slow on larger data sets.

- printAllBlocks:

  If `object` is a `BlockedInfinitySparseMatrix`, should summaries of
  all blocks be printed alongside the overall summary? Default `FALSE`.

- blockStructure:

  If `object` is a `BlockedInfinitySparseMatrix` and `printAllBlocks` is
  false, print a quick summary of each individual block. Default `TRUE`.
  If the number of blocks is high, consider suppressing this.

## Value

A named `list`. The summary for an `InfinitySparseMatrix` or
`DenseMatrix` contains the following:

- `total`: :

  Contains the total number of treatment and control members, as well as
  eligible and ineligible matches.

- `matchable`: :

  The names of all treatment and control members with at least one
  eligible match.

- `unmatchable`: :

  The names of all treatment and control members with no eligible
  matches.

- `distances`: :

  The summary of minimum matchable distances, if `distanceSummary` is
  `TRUE`.

For `BlockedInfinitySparseMatrix`, the named `list` instead of contains
one entry per block, named after each block (i.e. the value of the
blocking variable) as well as a block named 'overall' which contains the
summary ignoring blocks. Each of these entries contains a `list` with
entries 'total', 'matchable', 'unmatchable' and 'distances', as
described above.

## Details

The output consists of several pieces.

- Membership: Indicates the dimension of the distance.

- Total (in)eligible potential matches: A measure of the sparsity of the
  distance. Eligible matches have a finite distance between treatment
  and control members; they could be matched. Ineligible matches have
  `Inf` distance and can not be matched. A higher number of ineligible
  matches can speed up matching, but runs the risk of less optimal
  overall matching results.

- Unmatchable treatment/control members: If any observations have no
  eligible matches (e.g. their distance to every potential match is
  `Inf`) they are listed here. See Value below for details of how to
  access lists of matchable and unmatchable treatment and control
  members.

- Summary of minimum matchable distance per treatment member: To assist
  with choosing a caliper, this is a numeric summary of the smallest
  distance per matchable treatment member. If you provide a caliper that
  is less than the maximum value, at least one treatment member will
  become unmatchable.

- Block structure: For `BlockedInfinitySparseMatrix`, a quick summary of
  the structure of each individual block. (The above will all be across
  all blocks.) This may indicate which blocks, if any, are problematic.
