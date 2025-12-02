# Package index

## All functions

- [`BlockedInfinitySparseMatrix-class`](https://markmfredrickson.github.io/optmatch/dev/reference/BlockedInfinitySparseMatrix-class.md)
  : Blocked Infinity Sparse Matrix

- [`InfinitySparseMatrix-class`](https://markmfredrickson.github.io/optmatch/dev/reference/InfinitySparseMatrix-class.md)
  : Objects for sparse matching problems.

- [`LEMON()`](https://markmfredrickson.github.io/optmatch/dev/reference/LEMON.md)
  : (Internal) Helper function for accessing algorithms in LEMON solver

- [`antiExactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/antiExactMatch.md)
  : Specify a matching problem where units in a common factor cannot be
  matched.

- [`as.InfinitySparseMatrix()`](https://markmfredrickson.github.io/optmatch/dev/reference/as.InfinitySparseMatrix.md)
  : Convert an object to InfinitySparseMatrix

- [`as.list(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/as.list.BlockedInfinitySparseMatrix.md)
  : Splits a BlockedInfinitySparseMatrix into a list of
  InfinitySparseMatrices

- [`c(`*`<optmatch>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/c.optmatch.md)
  : Combine Optmatch objects

- [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
  : Prepare matching distances suitable for matching within calipers.

- [`cbind(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/cbindrbind.md)
  [`rbind(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/cbindrbind.md)
  [`cbind(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/cbindrbind.md)
  [`rbind(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/cbindrbind.md)
  : Combine InfinitySparseMatrices or BlockedInfinitySparseMatrices by
  row or column

- [`compare_optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/compare_optmatch.md)
  : Compares the equality of optmatch objects, ignoring attributes and
  group names.

- [`dbind()`](https://markmfredrickson.github.io/optmatch/dev/reference/dbind.md)
  : Diagonally bind together subgroup-specific distances

- [`dimnames(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/dimnames-InfinitySparseMatrix.md)
  [`` `dimnames<-`( ``*`<InfinitySparseMatrix>`*`,`*`<list>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/dimnames-InfinitySparseMatrix.md)
  [`` `dimnames<-`( ``*`<InfinitySparseMatrix>`*`,`*`<NULL>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/dimnames-InfinitySparseMatrix.md)
  : Get and set dimnames for InfinitySparseMatrix objects

- [`distUnion()`](https://markmfredrickson.github.io/optmatch/dev/reference/distUnion.md)
  : Combine multiple distance specifications into a single distance
  specification.

- [`effectiveSampleSize()`](https://markmfredrickson.github.io/optmatch/dev/reference/effectiveSampleSize.md)
  : Compute the effective sample size of a match.

- [`evaluate_primal()`](https://markmfredrickson.github.io/optmatch/dev/reference/evaluate_primal.md)
  : Compute value of primal problem given flows and arc costs

- [`exactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md)
  : Generate an exact matching set of subproblems.

- [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  : Create missingness indicator variables and non-informatively fill in
  missing values

- [`findSubproblems()`](https://markmfredrickson.github.io/optmatch/dev/reference/findSubproblems.md)
  : List subproblems of a distance

- [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  [`full()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  : Optimal full matching

- [`getMaxProblemSize()`](https://markmfredrickson.github.io/optmatch/dev/reference/getMaxProblemSize.md)
  : What is the maximum allowed problem size?

- [`subset(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ism.subset.md)
  [`` `[`( ``*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ism.subset.md)
  [`` `[<-`( ``*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ism.subset.md)
  : Subsetting for InfinitySparseMatrices

- [`` `+`( ``*`<InfinitySparseMatrix>`*`,`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ismBinaryOps.md)
  [`` `-`( ``*`<InfinitySparseMatrix>`*`,`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ismBinaryOps.md)
  [`` `*`( ``*`<InfinitySparseMatrix>`*`,`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ismBinaryOps.md)
  [`` `/`( ``*`<InfinitySparseMatrix>`*`,`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/ismBinaryOps.md)
  : Element-wise addition

- [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  : Create treated to control distances for matching problems

- [`matched()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md)
  [`unmatched()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md)
  [`matchfailed()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md)
  : Identification of units placed into matched sets

- [`matched.distances()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.distances.md)
  : Determine distances between matched units

- [`maxCaliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/maxCaliper.md)
  : Find the maximum caliper width that will create a feasible problem.

- [`c(`*`<SubProbInfo>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/mcf_c_fns.md)
  [`c(`*`<NodeInfo>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/mcf_c_fns.md)
  [`c(`*`<ArcInfo>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/mcf_c_fns.md)
  [`c(`*`<MCFSolutions>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/mcf_c_fns.md)
  [`c(`*`<FullmatchMCFSolutions>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/mcf_c_fns.md)
  : Combine objects

- [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  :

  (Deprecated, in favor of `match_on`) Create matching distances

- [`minExactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/minExactMatch.md)
  : Find the minimal exact match factors that will be feasible for a
  given maximum problem size.

- [`maxControlsCap()`](https://markmfredrickson.github.io/optmatch/dev/reference/minmaxctlcap.md)
  [`minControlsCap()`](https://markmfredrickson.github.io/optmatch/dev/reference/minmaxctlcap.md)
  : Set thinning and thickening caps for full matching

- [`nuclearplants`](https://markmfredrickson.github.io/optmatch/dev/reference/nuclearplants.md)
  : Nuclear Power Station Construction Data

- [`num_eligible_matches()`](https://markmfredrickson.github.io/optmatch/dev/reference/num_eligible_matches-methods.md)
  : Returns the number of eligible matches for the distance.

- [`pscore.dist()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch-defunct.md)
  [`mahal.dist()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch-defunct.md)
  : Functions deprecated or removed from optmatch

- [`summary(`*`<optmatch>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  : Optmatch Class

- [`optmatch_restrictions()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch_restrictions.md)
  : optmatch_restrictions

- [`optmatch_same_distance()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch_same_distance.md)
  :

  Checks if two distances are equivalent. `x` and `y` can be distances
  (`InfinitySparseMatrix`, `BlockedInfinitySparseMatrix`, or
  `DenseMatrix`), or they can be `optmatch` objects.

- [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  [`pair()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  : Optimal 1:1 and 1:k matching

- [`plantdist`](https://markmfredrickson.github.io/optmatch/dev/reference/plantdist.md)
  : Dissimilarities of Some U.S. Nuclear Plants

- [`predict(`*`<CBPS>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/predict.CBPS.md)
  : (Internal) Predict for CBPS objects

- [`print(`*`<optmatch>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/print.optmatch.md)
  :

  Printing `optmatch` objects.

- [`scoreCaliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/scoreCaliper.md)
  : (Internal) Helper function to create an InfinitySparseMatrix from a
  set of scores, a treatment indicator, and a caliper width.

- [`scores()`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
  : Extract scores (propensity, prognostic,...) from a fitted model

- [`setMaxProblemSize()`](https://markmfredrickson.github.io/optmatch/dev/reference/setMaxProblemSize.md)
  : Set the maximum problem size

- [`show(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/show-BlockedInfinitySparseMatrix-method.md)
  : Displays a BlockedInfinitySparseMatrix

- [`show(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/show-InfinitySparseMatrix-method.md)
  : Displays an InfinitySparseMatrix

- [`sort(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/sort.ism.md)
  [`sort(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/sort.ism.md)
  : Sort the internal structure of an InfinitySparseMatrix.

- [`strata()`](https://markmfredrickson.github.io/optmatch/dev/reference/strata.md)
  : Identify Stratafication Variables

- [`stratumStructure()`](https://markmfredrickson.github.io/optmatch/dev/reference/stratumStructure.md)
  [`print(`*`<stratumStructure>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/stratumStructure.md)
  : Return structure of matched sets

- [`subdim()`](https://markmfredrickson.github.io/optmatch/dev/reference/subdim-methods.md)
  : Returns the dimension of each valid subproblem

- [`summary(`*`<InfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md)
  [`summary(`*`<BlockedInfinitySparseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md)
  [`summary(`*`<DenseMatrix>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md)
  : Summarize a distance matrix

- [`update(`*`<optmatch>`*`)`](https://markmfredrickson.github.io/optmatch/dev/reference/update.optmatch.md)
  :

  Performs an update on an `optmatch` object.
