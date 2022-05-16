## Changes in **optmatch** Version 0.10.3

### Interface changes

- The rank Mahalanobis distance being created in `match_on()` using the argument
  `method = "rank_mahalanobis"` was accidentally returning the squared distance
  rather than the distance. This has been fixed. To recover results using the
  squared distance, square the results, e.g.: `match_on(..., method =
  "rank_mahalanobis")^2`. (Thanks Noah Greifer #218)
- New function `as.list.BlockedInfinitySparseMatrix()` to split a single
  `BlockedInfinitySparseMatrix` into a `list` of `InfinitySparseMatrix` based
  upon the separate blocks. (Called via `as.list(b)` when `b` is a
  `BlockedInfinitySparseMatrix`.)
- New function `dbind()` for binding several distance matrices into a single
  `BlockedInfinitySparseMatrix`. Valid inputs include any distance
  convertible into an `InfinitySparseMatrix`, or `BlockedInfinitySparseMatrix`,
  or `list`s of these. (#65)

### Infrastructure changes

- Manually correcting `License_is_FOSS` and `License_restricts_use` flags after
  0.10.0 transition to an open license.
- **optmatch** no longer depends on **digest**, **RItools**, or **survey**, or
  imports from **survival**. This should harden us against unexpected downtime
  should any of these packages be removed from CRAN.
  - Implemented `optmatch::strata` to be used in place of `survival::strata`.
    Loading **survival** and masking `strata` should not cause issues either.
  - Hashing of distance matrixes is now done internally.
  - **RItools** and **survery** are now only suggested with appropriate
      warnings if users attempt to utilize code with them without first
      installing the packages.
- Modernized some vignettes

(Note: 0.10.1 and 0.10.2 were functionally equivalent releases updated to
address an issue with CRAN and the `License_is_FOSS` and `License_restricts_use`
flags.)

## Changes in **optmatch** Version 0.10.0

### Major changes

- **optmatch** no longer includes the RELAX-IV solver internally. That solver
  can be still be used by installing the new **rrelaxiv** package (which will
  *not* be hosted on CRAN). When **rrelaxiv** is not available, **optmatch**
  instead uses a min-cost flow solver provided by the LEMON project,
  <https://lemon.cs.elte.hu/trac/lemon>; bindings to these are provided by the
  **rlemon** package. The user interface remains the same, other than a new
  optional argument for specifying which solver to use. In *very* limited
  testing, we've seen similar matching results with the new solver.
- The **rlemon** package offers bindings to four LEMON solvers. See
  `help(fullmatch)` for a discussion on those, and the new argument to
  `fullmatch()`, `solver = `.
- To continue using the RELAX-IV solver, whether for back-compatibility or as a
  matter of preference, install the **rrelaxiv** package, from
  <https://errickson.net/rrelaxiv/>. If **rrelaxiv** is installed, RELAX-IV will
  become the default solver automatically.

### Minor changes

- Remove dependence on the **digest** package when generating hashes of distance
  matrices.
- **RItools** is moved to Suggests instead of Imports.

## Changes in **optmatch** Version 0.9-17

- Fix to FORTRAN to conform with Writing R Extensions §6.6.2.
- Observations with NAs in a blocking variable are now retained, although marked
  as unmatchable (#206). And similarly for observations with NAs in a scalar
  matching variable (#189).
- Minor bug fix(es) incl. #211, #204

## Changes in **optmatch** Version 0.9-16
- Bug fix: integer overflow issue arising with large problems (#209)
- Minor refinement of support for glms from `survey::svyglm()` (#194)


## Changes in **optmatch** Version 0.9-15
- Small bug fix related to `survey::mad` and `survey::med` interface

## Changes in **optmatch** Version 0.9-14

- Bug fix: `within=` arguments to `match_on()`, or functions calling
  `match_on()` such as `pairmatch()` or `fullmatch()`, were sometimes ignored
  (#181).
- Binary operations on sparse matching distances now compare dimnames of the two
  before proceeding (#190).
- Bug fix: when `fullmatch()` or `pairmatch()` found it infeasible to create
  matches within an exact matching category, under some circumstances all
  members of that category were being placed into a single category labeled
  `1.NA`, or `2.NA` etc. Instead, all members of that category are now `NA`
  (#203).
- Fixed a bug causing `match_on()`, `scores()` to misinterpret propensity or
  other scores fitted with `survey::svyglm()` (#194).
- `boxplot()` gains a method for `svyglm` objects, e.g. propensity score models
  fitted with case weights via `survey::svyglm()` (#194).
- The meaning of one of `match_on.glm()`'s arguments has changed slightly: to
  circumvent scale standardization when matching on a propensity score or other
  index, you should now pass `standardization.scale = 1`, not
  `standardization.scale = NULL` (#194).

## Changes in **optmatch** Version 0.9-12

- Fixed a bug causing `summary.optmatch()` to fail b/c of NAs in the treatment
  variable (#155).
- Fixed bug with using custom distance functions (#180) and updated
  documentation related to custom distance functions.
- Fixed minor compatibility error in R-devel for math operations on sparse
  distance matrices (#179).

## Changes in **optmatch** Version 0.9-11
- Added `exclude` argument to `match_on()` mirroring the `exclude` argument for
  `caliper()`.
- `Optmatch` objects now support an `update()` function, `update.Optmatch()`.
  (#54)
- `Optmatch` objects can be combined via a `c()` function, `c.Optmatch()`. (#68)
- Added support for `labelled` treatment vectors which often arise when
  importing from Stata or SPSS. (#159)
- Introduced more informative error messages in a few situations. (#149, #104)
- Better handling of NA's in variables involved in
  matching/calipering/exactMatching. (#147)
- Fixed a bug with incorrect results in `matchfailed()`. (#175)

## Changes in **optmatch** Version 0.9-10
- Minor release to fix warnings during CRAN checks.

## Changes in **optmatch** Version 0.9-9
- Fixed a bug that caused the effective sample size to be rounded too
  aggresively in `summary.optmatch()`.
- Improved several error messages and warnings. (#138, #149, #142)
- Fixed use of `if(vectorOfThings)` usage that will give an error in upcoming R
  release.

## Changes in **optmatch** Version 0.9-8
- If pairmatch is asked to match within a stratum with fewer eligible controls
  than `controls` times the number of treatments, it now attempts to match in
  that stratum by leaving out some of the treatment units. (#116)
- The treatment indicator must be either numeric 0/1 (1 for treatment, 0 for
  control) or logical (TRUE for treatment, FALSE for control).
- Support for any other type of treatment vector (factors, character, other
  numerical values) has been deprecated. You can easily update treatment vectors
  using conditional statements, e.g. if you have a character "T" and "C",
  `treatment_new = treatment == "T"`.
- The treatment indicator can now include NA's. Any observations with NA
  treatment status will be excluded from distances matrices and will never
  match. WARNING: if the `data` argument is excluded from `fullmatch()` or
  `pairmatch()` and `num_NA` > 0 entries in the treatment status vector are NA,
  then the length of the vector produced by `fullmatch()` or `pairmatch()` won't
  match the length of the treatment status vector, having `num_NA` fewer
  observations. Don't forget to pass a `data` argument!
- Fixed bug affecting rank Mahalanobis matching in combination with calipers
  and/or exact matching constraints (#128)
- Addressed an issue affecting certain problems with both exact matching and
  caliper restrictions, where `min.controls`/`mean.controls`/`max.controls`
  directives would have been mistakenly applied to the wrong subclasses,
  resulting in strange warnings and, potentially, spurious match failures or
  unintended structural restrictions in some subclasses (#129).
- Structural restrictions allowing only many-one matches no longer cause
  `fullmatch()` to automatically fail. I.e. we've restored the behaviour of the
  software prior to version 0.8. (#132)

## Changes in **optmatch** Version 0.9-7
- Added support for **CBPS** created objects (#121).
- Improvments to documentation for several functions.
- Several small bugfixes.

## Changes in **optmatch** Version 0.9-6
- New material in vignettes, on general use of the package and on import/export
  of matching results and material between R and SAS or Stata (Josh Errickson).
- New `summary()` methods for `InfinitySparseMatrix`
  (`summary.InfinitySparseMatrix()`), `BlockedInfinitySparseMatrix`,
  (`summary.BlockedInfinitySparseMatrix()`) and `DenseMatrix`
  (`summary.DenseMatrix()`). I.e., you can call `summary()` on the result of a
  call to `match_on()` or `caliper()`. The information this returns may be
  useful for selecting caliper widths, and for managing computational burdens
  with large matching problems.
- Streamlined combinations of exact and propensity score matching. If you
  include "+ strata(fac)" on the right hand side of a propensity scoring model
  formula, then pass the fitted model to `pairmatch()`, `fullmatch()` or
  `match_on()`, then the factor "fac" will both serve as an independent variable
  for the propensity model and an exact matching variable (#101). See the
  examples on the help documentation for `fullmatch()`.
- `pairmatch()` and `fullmatch()` no longer generate "matched.distances"
   attributes for their results. To get this information, use
   `matched.distances()`.
- (Internal) methods for sorting of InfinitySparseMatrix's
- Deprecated: support passing the results of `fill.NAs()` directly to `glm()` or
  similar. Use the traditional formula and `data` argument version. See help
  documentation for `fill.NAs()` for examples.
- Fixed: Rcpp incompatibilities for some OSX users (4bbcaca); `boxplot()` method
  for fitted propensities ignoring varwidth argument (#113); various minor
  issues affecting package development and deployment (#110,...).

## Changes in **optmatch** Version 0.9-5
- Documentation adjustments.
- Explicit print method for output from explicit calls to `stratumStructure()`.

## Changes in **optmatch** Version 0.9-4

- Significant speed up of math operations for sparse distance objects (by Josh
  Buckner).
- Introducing `contr.match_on()`, a new default contrasts function for making
  Mahalanobis and Euclidean distances. Previously we used R defaults, which (a)
  generated different answers for the same factor depending on the ordering of
  the levels and (b) led to different distances for {0,1}-valued numeric
  variables and two level factors. (#80)
- match_on now takes strata as element of formula. Now users can write:
    match_on(z ~ x1 + x2 + strata(exactMatchVar)) Instead of match_on(z ~ x1 +
    x2, within = exactMatch(z ~ exactMatchVar))
- Fixed bug giving spurious infeasibility warnings, sample size reductions when
  using `fullmatch()` with feasible combinations of `min.controls`,
  `mean.controls`/`max.controls` and `max.controls` (#92)
- Various small bug fixes and documentation improvements.

## Changes in **optmatch** Version 0.9-3

- Fixed memory issues, potential segfaults in solver code. (Thank you, Peter
  Solenberger).
- Fixed bug in dropping cases with extraneous NAs when using `fullmatch()` or
  `pairmatch()` to create distance specifications directly.
- Fixed bug (#83) in `glm()` method for `match_on()` that caused observations
  with fixable NAs to be dropped too often.
- New function `distUnion()` allows combining arbitrary distance specifications.
- New function `antiExactMatch()` provides for matches that may only occur
  between treated and control units with _different_ values on a factor
  variable. This is the opposite of `exactMatch()`, which ensures matches occur
  within factor levels.
- Can now infer `data` argument in more cases when using the `summary()` method
  when the **RItools** package is present.
- Additional warnings and clarifications.

## Changes in **optmatch** Version 0.9-2

- Fixed issue #74 by properly setting the `omit.fraction` argument when there
  are unmatched controls.
- Improvements to the `minExactMatch()` function.
- Added `optmatch_verbose_message` option to provide additional warnings.
- Fixed crash when all NULL or NA vectors passed as arguments to `fullmatch()`.
- Added argument to `caliper()` function that allows returning values that fit
  the caliper instead of just indicators of which entries fit the caliper width.
- Calipers widths can be given per-treated unit, instead of globally.
- Additional binary operators for sparse matrix representations.
- Added new ranked Mahalanobis method for the formula method of `match_on()`.

## Changes in **optmatch** Version 0.9-1

- Subsetting of `Optmatch` objects now preserves (and subsets) the subproblem
  attribute.
- Performance improvements for match_on applied to glm's.
- The solver update of version 0.9-0 had a bug that in some circumstances caused
  hangups or malloc's [Issue #70]. We believe this is now fixed -- but please
  notify maintainer if you continue to experience the problem. (If you do, we'll
  reward you with an easy workaround.)

## Changes in **optmatch** Version 0.9-0

### NEW FEATURES

- Solver limits now depend on machine limits, not arbitrary constants defined by
  the **optmatch** maintainers. For large problems, users will see a warning,
  but the solver will attempt to solve.

- `fullmatch()` and `pairmatch()` can now take distance generating arguments
  directly, instead of having to first call `match_on()`. See the documentation
  for these two functions for more details.

- Infeasibility recovery in `fullmatch()`. When passing a combination of
  constraints (e.g. `max.controls`) that would make the matching infeasible,
  `fullmatch()` will now attempt to find a feasible match that respects those
  constraints, which will likely result in omitting some controls units.

- An additional argument to `fullmatch()`, `mean.controls`, is an alternative to
  the previous `omit.fraction`. (Only one of the two arguments can be
  presented.) The match will attempt to average mean.controls number of controls
  per treatment.

- Each `Optmatch` object now carries with it the constraints used to generate it
  (e.g. `max.controls`) as well as a hashed version of the distance it matched
  up, to help with some debugging/error checking but avoiding having to carry
  the entire distance matrix around.

- Creating a distance matrix prior to matching is now optional. `fullmatch()`
  now accepts arguments from which `match_on()` would create a distance, and
  create the match behind the scenes.

- Performance enhancements for distance calculations.

- Several new utility functions, including `subdim()`,
  `optmatch_restrictions()`, `optmatch_same_distance()`,
  `num_eligible_matches()`. See their help documentation for additional details.

- Arithmetic operations between InfinitySparseMatrices and vectors are
  supported. The operation is carried out as column by vector steps.

- `scores()` function allows including model predictions (such as propensity
  scores) in formulas directly (such as combining multiple propensity scores).
  The `scores()` function is preferred to predict() as it makes several smart
  choices to avoid dropping observations due to partial missingness and other
  useful preparations for matching.

### BUG FIXES

- `match_on()` is now a S3 generic function, which solves several bugs using
  propensity models from other packages.

- `summary()` method was giving overly pessimistic warnings about failures.

- fixed bug in how `Optmatch` objects were printing.

### DEPRECATED AND DEFUNCT
- `mdist()` is now deprecated, in favor of `match_on()`.

## Changes in **optmatch** Version 0.8-3
- Changes to make examples compatible with PDF manual

## Changes in **optmatch** Version 0.8-2

- `full()` and `pair()` are now aliases to `fullmatch()` and `pairmatch()`

- All `match_on()` methods take `caliper` arguments (formerly just the numeric
  method and derived methods had this argument).

- boxplot methods for fitted propensity score methods (`glm()` and `bigglm()`)

- `fill.NAs()` now takes `contrasts.arg` argument to mimic `model.matrix()`

- Several bug fixes in examples, documentation

- The methods `pscore.dist()` and `mahal.dist()` are now deprecated, with useful
  error messages pointing users to replacements.

- Significant performance improvements for sparse matching problems.

- Functions `umatched()` and `matched()` were backwards. Corrected.

## Changes in **optmatch** Version 0.8-1

- Several small bug fixes

## Changes in **optmatch** Version 0.8-0

### NEW FEATURES

- More efficient data structure for sparse matching problems, those with
  relatively few allowed (finite) distances between units. Sparse problems often
  arise when calipers are employed. The new data structure
  (`InfinitySparseMatrix`) behaves like a simple matrix, allowing `cbind()`,
  `rbind()`, and `subset()` operations, making it easier to work with the older
  `optmatch.dlist` data structure.

- `match_on()`: A series of methods to generate matching problems using the new
  data structure when appropriate, or using a standard matrix when the problem
  is dense. This function is being deployed along side the `mdist()` function to
  provide complete backward compatibility. New development will focus on this
  function for distance creation, and users are encouraged to use it right away.
  One difference for `mdist()` users is the `within` argument. This argument
  takes an existing distance specification and limits the new comparisons to
  only those pairs that have finite distances in the `within` argument. See the
  `match_on()`, `exactMatch()`, and `caliper()` documentation for more details.

- `exactMatch()`: A new function to create stratified matching problems (in
  which cross strata matches are forbidden). Users can specify the strata using
  either a factor vector or a convenient formula interface. The results can be
  used in calls `match_on()` to limit distance calculations to only with-in
  strata treatment-control pairs.

- New `data` argument to `fullmatch()` and `pairmatch()`: This argument will set
  the order of the match to that of the `row.names`, `names`, or contents of the
  passed `data.frame` or `vector`. This avoids potential bugs caused when the
  `optmatch` objects were in a different order than users' data.

- Test suite expanded and now uses the **testthat** library.

- `fill.NAs()` allows (optionally) filling in all columns (previously, the first
  column was assumed to be an outcome or treatment indicator and was not filled
  in).

- New tools to find minimum feasible constraints: Large matching problems could
  exceed the upper limit for a matching problem. The functions `minExactmatch()`
  and `maxCaliper()` find the smallest interaction of potential factors for
  stratified matchings or the largest (most generous) caliper, respectively,
  that make the problem small enough to fit under the maximum problem size
  limit. See the help pages for these functions for more information.

### BUG FIXES

- Unmatched units are always NA (instead of being labeled `1.NA` or similar).
  This avoids some obscure bugs when feeding the results of `fullmatch()` to
  other functions.

FOR A DETAILED CHANGELOG, SEE <https://github.com/markmfredrickson/optmatch>

## Changes in **optmatch** Version 0.7-1

### NEW FEATURES

- `pairmatch()` has a new option, `remove.unmatchables`, that may be useful in
  conjunction with caliper matching. With `remove.unmatchables = TRUE`, prior to
  matching any units with no counterparts within caliper distance are removed.
  Pair matching can still fail, if for example for two distinct treatment units
  only a single control, the same one, is available for matching to them; but
  `remove.unmatchables` eliminates one simple and common reason for pair
  matching to fail.

- Applying `summary()` to an optmatch object now creates a `summary.optmatch`
  containing the summary information, in addition to reporting it to the console
  (via a `summary.optmatch()` method for `print()`).

- `mdist.formula()` no longer requires an explicit data argument. I.e., you can
  get away with a call like `mdist(Treat~X1+X2|S)` if the variables `Treat`,
  `X1`, `X2` and `S` are available in the environment you're working from (or in
  one of its parent environments). Previously you would have had to do
  `mdist(Treat~X1+X2|S, data=mydata)`. (The latter formulation is still to be
  preferred, however, in part because with it `mdist()` gets to use data's row
  names, whereas otherwise it would have to make up row names.)

## Changes in **optmatch** Version 0.7

### NEW FEATURES

- New function `fill.NAs()` replaces missing observations (ie. NA values) with
  minimally informative values (ie. the mean of observed columns). `fill.NAs()`
  handles functions in formulas intelligently and provides missing indicators
  for each variable. See the help documentation for more information and
  examples.

### BUG FIXES

- `mdist.function()` method now properly returns an `optmatch.dlist` object for
  use in `summary.optmatch()`, etc.

- `mdist.function()` maintains label on grouping factor.

## Changes in **optmatch** Version 0.6-1

### NEW FEATURES

- New `mdist()` method to extract propensity scores from models fitted using
  `bigglm()` in package **biglm**.

- `mdist()`'s formula method now understands grouping factors indicated with a
  pipe (`|`)

- informative error message for `mdist()` called on numeric vectors

- updated `mdist()` documentation

## Changes in **optmatch** Version 0.6

### NEW FEATURES

- There is a new generic function, `mdist()`, for creating matching distances.
  It accepts: fitted glm's, which it uses to extract propensity distances;
  formulas, which it uses to construct squared Mahalanobis distances; and
  functions, with which a user can construct his or her own type of distance.
  The function method is more intuitive to work with than the older `makedist()`
  function.

- A new function, `caliper()`, builds on the `mdist()` structure to provide a
  convenient way to add calipers to a distance. In contrast to earlier ways of
  adding calipers, `caliper()` has an optional argument specify observations to
  be excluded from the caliper requirement --- this permits one to relax it for
  just a few observations, for instance.

- `summary.optmatch()` now removes strata in which matching failed (b/c the
  matching problem was found to be infeasible) before summarizing. It also
  indicates when such strata are present, and how many observations fall in
  them.

- Demo has been updated to reflect changes as of version 0.4, 0.5, 0.6.

### DEPRECATED & DEFUNCT

- The vignette is sufficiently out of date that it's been removed.

### BUG FIXES

- subsetting of objects of class `Optmatch` now preserves matched.distances
  attribute.

- fixed bug in `maxControlsCap()`/`minControlsCap()` whereby they behaved
  unreliably on subclasses within which some subjects had no permissible
  matches.

- Removed unnecessary panic in `fullmatch()` when it was given a `min.controls`
  argument with attributes other than names (as when it is created by
  `tapply()`).

- fixed bug wherein `summary.optmatch()` fails to retrieve balance tests if
  given a propensity model that had function calls in its formula.

- Documentation pages for `fullmatch()`, `pairmatch()` filled out a bit.

## Changes in **optmatch** Version 0.5

### NEW FEATURES:

- `summary.optmatch()` completely revised. It now reports information about the
  configuration of the matched sets and about matched distances. In addition, if
  given a fitted propensity model as a second argument it summarizes covariate
  balance.
