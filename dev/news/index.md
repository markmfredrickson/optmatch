# Changelog

## Changes in **optmatch** Version 0.10.8

CRAN release: 2024-09-19

- `optmatch:::scoreCaliper()` gains an optional `within=` argument
  ([\#245](https://github.com/markmfredrickson/optmatch/issues/245))
- Updates to internal C++ code

## Changes in **optmatch** Version 0.10.7

CRAN release: 2023-11-15

- Hardened tests further against package unavailability
- Changed the name of default arguments for a version of `predict`.
  ([\#223](https://github.com/markmfredrickson/optmatch/issues/223))
- Several small documentation tweaks to pass CRAN checks on R-devel.
- Fix a bug introduced in version 0.10.6, involving discretization of
  distances. The fix avoids spurious errors for distance matrices with
  very large values, although you may still have to pass tol= arguments
  to pairmatch() and fullmatch() that are smaller than the desired
  tolerance.([\#230](https://github.com/markmfredrickson/optmatch/issues/230))
- Disable passing of local variables from generic corresponding to
  change in R_USEMETHOD_FORWARD_LOCALS coming in the next major release.
  No user-facing change, except the `@call` slot of objects may look
  slighly different (but should function identically).
  ([\#234](https://github.com/markmfredrickson/optmatch/issues/234))

## Changes in **optmatch** Version 0.10.6

CRAN release: 2023-02-07

- Adjusted a check to avoid ambiguous failures when using LEMON vs
  RELAX-IV.
- Updated CITATION to use
  [`bibentry()`](https://rdrr.io/r/utils/bibentry.html).
- Minor tweaks to address failing tests.

## Changes in **optmatch** Version 0.10.5

CRAN release: 2022-08-15

Minor revision to address a failing test.

## Changes in **optmatch** Version 0.10.4

CRAN release: 2022-08-09

- When including factor variables on the right hand side of the formula
  passed into
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  now more simply calculates the contrast to enable more intuitive
  results. (Thanks Noah Greifer,
  [\#220](https://github.com/markmfredrickson/optmatch/issues/220))
- [`dbind()`](https://markmfredrickson.github.io/optmatch/dev/reference/dbind.md)
  will now properly support binding more than 26 unique matrices when
  renaming is necessary; in fact it supports up to 18,278 uniquely
  renamed matrices.
- A few tweaks in documentation, testing and vignettes to satisfy CRAN
  requirements and harden against loss of dependencies.

## Changes in **optmatch** Version 0.10.3

CRAN release: 2022-05-16

### Interface changes

- The rank Mahalanobis distance being created in
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  using the argument `method = "rank_mahalanobis"` was accidentally
  returning the squared distance rather than the distance. This has been
  fixed. To recover results using the squared distance, square the
  results, e.g.: `match_on(..., method = "rank_mahalanobis")^2`. (Thanks
  Noah Greifer
  [\#218](https://github.com/markmfredrickson/optmatch/issues/218))
- New function
  [`as.list.BlockedInfinitySparseMatrix()`](https://markmfredrickson.github.io/optmatch/dev/reference/as.list.BlockedInfinitySparseMatrix.md)
  to split a single `BlockedInfinitySparseMatrix` into a `list` of
  `InfinitySparseMatrix` based upon the separate blocks. (Called via
  `as.list(b)` when `b` is a `BlockedInfinitySparseMatrix`.)
- New function
  [`dbind()`](https://markmfredrickson.github.io/optmatch/dev/reference/dbind.md)
  for binding several distance matrices into a single
  `BlockedInfinitySparseMatrix`. Valid inputs include any distance
  convertible into an `InfinitySparseMatrix`, or
  `BlockedInfinitySparseMatrix`, or `list`s of these.
  ([\#65](https://github.com/markmfredrickson/optmatch/issues/65))

### Infrastructure changes

- Manually correcting `License_is_FOSS` and `License_restricts_use`
  flags after 0.10.0 transition to an open license.
- **optmatch** no longer depends on **digest**, **RItools**, or
  **survey**, or imports from **survival**. This should harden us
  against unexpected downtime should any of these packages be removed
  from CRAN.
  - Implemented
    [`optmatch::strata`](https://markmfredrickson.github.io/optmatch/dev/reference/strata.md)
    to be used in place of
    [`survival::strata`](https://rdrr.io/pkg/survival/man/strata.html).
    Loading **survival** and masking `strata` should not cause issues
    either.
  - Hashing of distance matrixes is now done internally.
  - **RItools** and **survey** are now only suggested with appropriate
    warnings if users attempt to utilize code with them without first
    installing the packages.
- Modernized some vignettes

(Note: 0.10.1 and 0.10.2 were functionally equivalent releases updated
to address an issue with CRAN and the `License_is_FOSS` and
`License_restricts_use` flags.)

## Changes in **optmatch** Version 0.10.0

CRAN release: 2022-03-15

### Major changes

- **optmatch** no longer includes the RELAX-IV solver internally. That
  solver can be still be used by installing the new **rrelaxiv** package
  (which will *not* be hosted on CRAN). When **rrelaxiv** is not
  available, **optmatch** instead uses a min-cost flow solver provided
  by the LEMON project, <https://lemon.cs.elte.hu/trac/lemon>; bindings
  to these are provided by the **rlemon** package. The user interface
  remains the same, other than a new optional argument for specifying
  which solver to use. In *very* limited testing, we’ve seen similar
  matching results with the new solver.
- The **rlemon** package offers bindings to four LEMON solvers. See
  [`help(fullmatch)`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  for a discussion on those, and the new argument to
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
  `solver =`.
- To continue using the RELAX-IV solver, whether for back-compatibility
  or as a matter of preference, install the **rrelaxiv** package, from
  <https://errickson.net/rrelaxiv/>. If **rrelaxiv** is installed,
  RELAX-IV will become the default solver automatically.

### Minor changes

- Remove dependence on the **digest** package when generating hashes of
  distance matrices.
- **RItools** is moved to Suggests instead of Imports.

## Changes in **optmatch** Version 0.9-17

CRAN release: 2022-02-22

- Fix to FORTRAN to conform with Writing R Extensions §6.6.2.
- Observations with NAs in a blocking variable are now retained,
  although marked as unmatchable
  ([\#206](https://github.com/markmfredrickson/optmatch/issues/206)).
  And similarly for observations with NAs in a scalar matching variable
  ([\#189](https://github.com/markmfredrickson/optmatch/issues/189)).
- Minor bug fix(es) incl. #211,
  [\#204](https://github.com/markmfredrickson/optmatch/issues/204)

## Changes in **optmatch** Version 0.9-16

- Bug fix: integer overflow issue arising with large problems
  ([\#209](https://github.com/markmfredrickson/optmatch/issues/209))
- Minor refinement of support for glms from
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  ([\#194](https://github.com/markmfredrickson/optmatch/issues/194))

## Changes in **optmatch** Version 0.9-15

CRAN release: 2021-08-17

- Small bug fix related to `survey::mad` and `survey::med` interface

## Changes in **optmatch** Version 0.9-14

CRAN release: 2021-05-26

- Bug fix: `within=` arguments to
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  or functions calling
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  such as
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  or
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
  were sometimes ignored
  ([\#181](https://github.com/markmfredrickson/optmatch/issues/181)).
- Binary operations on sparse matching distances now compare dimnames of
  the two before proceeding
  ([\#190](https://github.com/markmfredrickson/optmatch/issues/190)).
- Bug fix: when
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  found it infeasible to create matches within an exact matching
  category, under some circumstances all members of that category were
  being placed into a single category labeled `1.NA`, or `2.NA` etc.
  Instead, all members of that category are now `NA`
  ([\#203](https://github.com/markmfredrickson/optmatch/issues/203)).
- Fixed a bug causing
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  [`scores()`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
  to misinterpret propensity or other scores fitted with
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  ([\#194](https://github.com/markmfredrickson/optmatch/issues/194)).
- [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) gains a method
  for `svyglm` objects, e.g. propensity score models fitted with case
  weights via
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  ([\#194](https://github.com/markmfredrickson/optmatch/issues/194)).
- The meaning of one of
  [`match_on.glm()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)’s
  arguments has changed slightly: to circumvent scale standardization
  when matching on a propensity score or other index, you should now
  pass `standardization.scale = 1`, not `standardization.scale = NULL`
  ([\#194](https://github.com/markmfredrickson/optmatch/issues/194)).

## Changes in **optmatch** Version 0.9-12

CRAN release: 2019-10-11

- Fixed a bug causing
  [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  to fail b/c of NAs in the treatment variable
  ([\#155](https://github.com/markmfredrickson/optmatch/issues/155)).
- Fixed bug with using custom distance functions
  ([\#180](https://github.com/markmfredrickson/optmatch/issues/180)) and
  updated documentation related to custom distance functions.
- Fixed minor compatibility error in R-devel for math operations on
  sparse distance matrices
  ([\#179](https://github.com/markmfredrickson/optmatch/issues/179)).

## Changes in **optmatch** Version 0.9-11

CRAN release: 2019-07-12

- Added `exclude` argument to
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  mirroring the `exclude` argument for
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md).
- `Optmatch` objects now support an
  [`update()`](https://rdrr.io/r/stats/update.html) function,
  `update.Optmatch()`.
  ([\#54](https://github.com/markmfredrickson/optmatch/issues/54))
- `Optmatch` objects can be combined via a
  [`c()`](https://rdrr.io/r/base/c.html) function, `c.Optmatch()`.
  ([\#68](https://github.com/markmfredrickson/optmatch/issues/68))
- Added support for `labelled` treatment vectors which often arise when
  importing from Stata or SPSS.
  ([\#159](https://github.com/markmfredrickson/optmatch/issues/159))
- Introduced more informative error messages in a few situations.
  ([\#149](https://github.com/markmfredrickson/optmatch/issues/149),
  [\#104](https://github.com/markmfredrickson/optmatch/issues/104))
- Better handling of NA’s in variables involved in
  matching/calipering/exactMatching.
  ([\#147](https://github.com/markmfredrickson/optmatch/issues/147))
- Fixed a bug with incorrect results in
  [`matchfailed()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md).
  ([\#175](https://github.com/markmfredrickson/optmatch/issues/175))

## Changes in **optmatch** Version 0.9-10

CRAN release: 2018-07-12

- Minor release to fix warnings during CRAN checks.

## Changes in **optmatch** Version 0.9-9

CRAN release: 2018-05-14

- Fixed a bug that caused the effective sample size to be rounded too
  aggresively in
  [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md).
- Improved several error messages and warnings.
  ([\#138](https://github.com/markmfredrickson/optmatch/issues/138),
  [\#149](https://github.com/markmfredrickson/optmatch/issues/149),
  [\#142](https://github.com/markmfredrickson/optmatch/issues/142))
- Fixed use of `if(vectorOfThings)` usage that will give an error in
  upcoming R release.

## Changes in **optmatch** Version 0.9-8

CRAN release: 2018-01-16

- If pairmatch is asked to match within a stratum with fewer eligible
  controls than `controls` times the number of treatments, it now
  attempts to match in that stratum by leaving out some of the treatment
  units.
  ([\#116](https://github.com/markmfredrickson/optmatch/issues/116))
- The treatment indicator must be either numeric 0/1 (1 for treatment, 0
  for control) or logical (TRUE for treatment, FALSE for control).
- Support for any other type of treatment vector (factors, character,
  other numerical values) has been deprecated. You can easily update
  treatment vectors using conditional statements, e.g. if you have a
  character “T” and “C”, `treatment_new = treatment == "T"`.
- The treatment indicator can now include NA’s. Any observations with NA
  treatment status will be excluded from distances matrices and will
  never match. WARNING: if the `data` argument is excluded from
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  and `num_NA` \> 0 entries in the treatment status vector are NA, then
  the length of the vector produced by
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  won’t match the length of the treatment status vector, having `num_NA`
  fewer observations. Don’t forget to pass a `data` argument!
- Fixed bug affecting rank Mahalanobis matching in combination with
  calipers and/or exact matching constraints
  ([\#128](https://github.com/markmfredrickson/optmatch/issues/128))
- Addressed an issue affecting certain problems with both exact matching
  and caliper restrictions, where
  `min.controls`/`mean.controls`/`max.controls` directives would have
  been mistakenly applied to the wrong subclasses, resulting in strange
  warnings and, potentially, spurious match failures or unintended
  structural restrictions in some subclasses
  ([\#129](https://github.com/markmfredrickson/optmatch/issues/129)).
- Structural restrictions allowing only many-one matches no longer cause
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  to automatically fail. I.e. we’ve restored the behaviour of the
  software prior to version 0.8.
  ([\#132](https://github.com/markmfredrickson/optmatch/issues/132))

## Changes in **optmatch** Version 0.9-7

CRAN release: 2016-12-29

- Added support for **CBPS** created objects
  ([\#121](https://github.com/markmfredrickson/optmatch/issues/121)).
- Improvments to documentation for several functions.
- Several small bugfixes.

## Changes in **optmatch** Version 0.9-6

CRAN release: 2016-05-03

- New material in vignettes, on general use of the package and on
  import/export of matching results and material between R and SAS or
  Stata (Josh Errickson).
- New [`summary()`](https://rdrr.io/r/base/summary.html) methods for
  `InfinitySparseMatrix`
  ([`summary.InfinitySparseMatrix()`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md)),
  `BlockedInfinitySparseMatrix`,
  ([`summary.BlockedInfinitySparseMatrix()`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md))
  and `DenseMatrix`
  ([`summary.DenseMatrix()`](https://markmfredrickson.github.io/optmatch/dev/reference/summary.ism.md)).
  I.e., you can call [`summary()`](https://rdrr.io/r/base/summary.html)
  on the result of a call to
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  or
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md).
  The information this returns may be useful for selecting caliper
  widths, and for managing computational burdens with large matching
  problems.
- Streamlined combinations of exact and propensity score matching. If
  you include “+ strata(fac)” on the right hand side of a propensity
  scoring model formula, then pass the fitted model to
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md),
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  then the factor “fac” will both serve as an independent variable for
  the propensity model and an exact matching variable
  ([\#101](https://github.com/markmfredrickson/optmatch/issues/101)).
  See the examples on the help documentation for
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md).
- [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  and
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  no longer generate “matched.distances” attributes for their results.
  To get this information, use
  [`matched.distances()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.distances.md).
- (Internal) methods for sorting of InfinitySparseMatrix’s
- Deprecated: support passing the results of
  [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  directly to [`glm()`](https://rdrr.io/r/stats/glm.html) or similar.
  Use the traditional formula and `data` argument version. See help
  documentation for
  [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  for examples.
- Fixed: Rcpp incompatibilities for some OSX users (4bbcaca);
  [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) method for
  fitted propensities ignoring varwidth argument
  ([\#113](https://github.com/markmfredrickson/optmatch/issues/113));
  various minor issues affecting package development and deployment
  ([\#110](https://github.com/markmfredrickson/optmatch/issues/110),…).

## Changes in **optmatch** Version 0.9-5

CRAN release: 2015-07-31

- Documentation adjustments.
- Explicit print method for output from explicit calls to
  [`stratumStructure()`](https://markmfredrickson.github.io/optmatch/dev/reference/stratumStructure.md).

## Changes in **optmatch** Version 0.9-4

- Significant speed up of math operations for sparse distance objects
  (by Josh Buckner).
- Introducing `contr.match_on()`, a new default contrasts function for
  making Mahalanobis and Euclidean distances. Previously we used R
  defaults, which (a) generated different answers for the same factor
  depending on the ordering of the levels and (b) led to different
  distances for {0,1}-valued numeric variables and two level factors.
  ([\#80](https://github.com/markmfredrickson/optmatch/issues/80))
- match_on now takes strata as element of formula. Now users can write:
  match_on(z ~ x1 + x2 + strata(exactMatchVar)) Instead of match_on(z ~
  x1 + x2, within = exactMatch(z ~ exactMatchVar))
- Fixed bug giving spurious infeasibility warnings, sample size
  reductions when using
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  with feasible combinations of `min.controls`,
  `mean.controls`/`max.controls` and `max.controls`
  ([\#92](https://github.com/markmfredrickson/optmatch/issues/92))
- Various small bug fixes and documentation improvements.

## Changes in **optmatch** Version 0.9-3

CRAN release: 2014-08-13

- Fixed memory issues, potential segfaults in solver code. (Thank you,
  Peter Solenberger).
- Fixed bug in dropping cases with extraneous NAs when using
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  to create distance specifications directly.
- Fixed bug
  ([\#83](https://github.com/markmfredrickson/optmatch/issues/83)) in
  [`glm()`](https://rdrr.io/r/stats/glm.html) method for
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  that caused observations with fixable NAs to be dropped too often.
- New function
  [`distUnion()`](https://markmfredrickson.github.io/optmatch/dev/reference/distUnion.md)
  allows combining arbitrary distance specifications.
- New function
  [`antiExactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/antiExactMatch.md)
  provides for matches that may only occur between treated and control
  units with *different* values on a factor variable. This is the
  opposite of
  [`exactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
  which ensures matches occur within factor levels.
- Can now infer `data` argument in more cases when using the
  [`summary()`](https://rdrr.io/r/base/summary.html) method when the
  **RItools** package is present.
- Additional warnings and clarifications.

## Changes in **optmatch** Version 0.9-2

- Fixed issue
  [\#74](https://github.com/markmfredrickson/optmatch/issues/74) by
  properly setting the `omit.fraction` argument when there are unmatched
  controls.
- Improvements to the
  [`minExactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/minExactMatch.md)
  function.
- Added `optmatch_verbose_message` option to provide additional
  warnings.
- Fixed crash when all NULL or NA vectors passed as arguments to
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md).
- Added argument to
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
  function that allows returning values that fit the caliper instead of
  just indicators of which entries fit the caliper width.
- Calipers widths can be given per-treated unit, instead of globally.
- Additional binary operators for sparse matrix representations.
- Added new ranked Mahalanobis method for the formula method of
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md).

## Changes in **optmatch** Version 0.9-1

CRAN release: 2013-11-01

- Subsetting of `Optmatch` objects now preserves (and subsets) the
  subproblem attribute.
- Performance improvements for match_on applied to glm’s.
- The solver update of version 0.9-0 had a bug that in some
  circumstances caused hangups or malloc’s \[Issue
  [\#70](https://github.com/markmfredrickson/optmatch/issues/70)\]. We
  believe this is now fixed – but please notify maintainer if you
  continue to experience the problem. (If you do, we’ll reward you with
  an easy workaround.)

## Changes in **optmatch** Version 0.9-0

CRAN release: 2013-08-10

### NEW FEATURES

- Solver limits now depend on machine limits, not arbitrary constants
  defined by the **optmatch** maintainers. For large problems, users
  will see a warning, but the solver will attempt to solve.

- [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  and
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  can now take distance generating arguments directly, instead of having
  to first call
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md).
  See the documentation for these two functions for more details.

- Infeasibility recovery in
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md).
  When passing a combination of constraints (e.g. `max.controls`) that
  would make the matching infeasible,
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  will now attempt to find a feasible match that respects those
  constraints, which will likely result in omitting some controls units.

- An additional argument to
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
  `mean.controls`, is an alternative to the previous `omit.fraction`.
  (Only one of the two arguments can be presented.) The match will
  attempt to average mean.controls number of controls per treatment.

- Each `Optmatch` object now carries with it the constraints used to
  generate it (e.g. `max.controls`) as well as a hashed version of the
  distance it matched up, to help with some debugging/error checking but
  avoiding having to carry the entire distance matrix around.

- Creating a distance matrix prior to matching is now optional.
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  now accepts arguments from which
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  would create a distance, and create the match behind the scenes.

- Performance enhancements for distance calculations.

- Several new utility functions, including
  [`subdim()`](https://markmfredrickson.github.io/optmatch/dev/reference/subdim-methods.md),
  [`optmatch_restrictions()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch_restrictions.md),
  [`optmatch_same_distance()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch_same_distance.md),
  [`num_eligible_matches()`](https://markmfredrickson.github.io/optmatch/dev/reference/num_eligible_matches-methods.md).
  See their help documentation for additional details.

- Arithmetic operations between InfinitySparseMatrices and vectors are
  supported. The operation is carried out as column by vector steps.

- [`scores()`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
  function allows including model predictions (such as propensity
  scores) in formulas directly (such as combining multiple propensity
  scores). The
  [`scores()`](https://markmfredrickson.github.io/optmatch/dev/reference/scores.md)
  function is preferred to predict() as it makes several smart choices
  to avoid dropping observations due to partial missingness and other
  useful preparations for matching.

### BUG FIXES

- [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  is now a S3 generic function, which solves several bugs using
  propensity models from other packages.

- [`summary()`](https://rdrr.io/r/base/summary.html) method was giving
  overly pessimistic warnings about failures.

- fixed bug in how `Optmatch` objects were printing.

### DEPRECATED AND DEFUNCT

- [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  is now deprecated, in favor of
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md).

## Changes in **optmatch** Version 0.8-3

CRAN release: 2013-05-16

- Changes to make examples compatible with PDF manual

## Changes in **optmatch** Version 0.8-2

- [`full()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  and
  [`pair()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  are now aliases to
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  and
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)

- All
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  methods take `caliper` arguments (formerly just the numeric method and
  derived methods had this argument).

- boxplot methods for fitted propensity score methods
  ([`glm()`](https://rdrr.io/r/stats/glm.html) and `bigglm()`)

- [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  now takes `contrasts.arg` argument to mimic
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)

- Several bug fixes in examples, documentation

- The methods
  [`pscore.dist()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch-defunct.md)
  and
  [`mahal.dist()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch-defunct.md)
  are now deprecated, with useful error messages pointing users to
  replacements.

- Significant performance improvements for sparse matching problems.

- Functions `umatched()` and
  [`matched()`](https://markmfredrickson.github.io/optmatch/dev/reference/matched.md)
  were backwards. Corrected.

## Changes in **optmatch** Version 0.8-1

CRAN release: 2013-01-03

- Several small bug fixes

## Changes in **optmatch** Version 0.8-0

CRAN release: 2012-11-17

### NEW FEATURES

- More efficient data structure for sparse matching problems, those with
  relatively few allowed (finite) distances between units. Sparse
  problems often arise when calipers are employed. The new data
  structure (`InfinitySparseMatrix`) behaves like a simple matrix,
  allowing [`cbind()`](https://rdrr.io/r/base/cbind.html),
  [`rbind()`](https://rdrr.io/r/base/cbind.html), and
  [`subset()`](https://rdrr.io/r/base/subset.html) operations, making it
  easier to work with the older `optmatch.dlist` data structure.

- [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md):
  A series of methods to generate matching problems using the new data
  structure when appropriate, or using a standard matrix when the
  problem is dense. This function is being deployed along side the
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  function to provide complete backward compatibility. New development
  will focus on this function for distance creation, and users are
  encouraged to use it right away. One difference for
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  users is the `within` argument. This argument takes an existing
  distance specification and limits the new comparisons to only those
  pairs that have finite distances in the `within` argument. See the
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md),
  [`exactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md),
  and
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
  documentation for more details.

- [`exactMatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/exactMatch.md):
  A new function to create stratified matching problems (in which cross
  strata matches are forbidden). Users can specify the strata using
  either a factor vector or a convenient formula interface. The results
  can be used in calls
  [`match_on()`](https://markmfredrickson.github.io/optmatch/dev/reference/match_on-methods.md)
  to limit distance calculations to only with-in strata
  treatment-control pairs.

- New `data` argument to
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  and
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md):
  This argument will set the order of the match to that of the
  `row.names`, `names`, or contents of the passed `data.frame` or
  `vector`. This avoids potential bugs caused when the `optmatch`
  objects were in a different order than users’ data.

- Test suite expanded and now uses the **testthat** library.

- [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  allows (optionally) filling in all columns (previously, the first
  column was assumed to be an outcome or treatment indicator and was not
  filled in).

- New tools to find minimum feasible constraints: Large matching
  problems could exceed the upper limit for a matching problem. The
  functions `minExactmatch()` and
  [`maxCaliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/maxCaliper.md)
  find the smallest interaction of potential factors for stratified
  matchings or the largest (most generous) caliper, respectively, that
  make the problem small enough to fit under the maximum problem size
  limit. See the help pages for these functions for more information.

### BUG FIXES

- Unmatched units are always NA (instead of being labeled `1.NA` or
  similar). This avoids some obscure bugs when feeding the results of
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  to other functions.

FOR A DETAILED CHANGELOG, SEE
<https://github.com/markmfredrickson/optmatch>

## Changes in **optmatch** Version 0.7-1

CRAN release: 2011-03-06

### NEW FEATURES

- [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  has a new option, `remove.unmatchables`, that may be useful in
  conjunction with caliper matching. With `remove.unmatchables = TRUE`,
  prior to matching any units with no counterparts within caliper
  distance are removed. Pair matching can still fail, if for example for
  two distinct treatment units only a single control, the same one, is
  available for matching to them; but `remove.unmatchables` eliminates
  one simple and common reason for pair matching to fail.

- Applying [`summary()`](https://rdrr.io/r/base/summary.html) to an
  optmatch object now creates a `summary.optmatch` containing the
  summary information, in addition to reporting it to the console (via a
  [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  method for [`print()`](https://rdrr.io/r/base/print.html)).

- [`mdist.formula()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  no longer requires an explicit data argument. I.e., you can get away
  with a call like `mdist(Treat~X1+X2|S)` if the variables `Treat`,
  `X1`, `X2` and `S` are available in the environment you’re working
  from (or in one of its parent environments). Previously you would have
  had to do `mdist(Treat~X1+X2|S, data=mydata)`. (The latter formulation
  is still to be preferred, however, in part because with it
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  gets to use data’s row names, whereas otherwise it would have to make
  up row names.)

## Changes in **optmatch** Version 0.7

### NEW FEATURES

- New function
  [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  replaces missing observations (ie. NA values) with minimally
  informative values (ie. the mean of observed columns).
  [`fill.NAs()`](https://markmfredrickson.github.io/optmatch/dev/reference/fill.NAs.md)
  handles functions in formulas intelligently and provides missing
  indicators for each variable. See the help documentation for more
  information and examples.

### BUG FIXES

- [`mdist.function()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  method now properly returns an `optmatch.dlist` object for use in
  [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md),
  etc.

- [`mdist.function()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  maintains label on grouping factor.

## Changes in **optmatch** Version 0.6-1

CRAN release: 2010-02-11

### NEW FEATURES

- New
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  method to extract propensity scores from models fitted using
  `bigglm()` in package **biglm**.

- [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)’s
  formula method now understands grouping factors indicated with a pipe
  (`|`)

- informative error message for
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  called on numeric vectors

- updated
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  documentation

## Changes in **optmatch** Version 0.6

### NEW FEATURES

- There is a new generic function,
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md),
  for creating matching distances. It accepts: fitted glm’s, which it
  uses to extract propensity distances; formulas, which it uses to
  construct squared Mahalanobis distances; and functions, with which a
  user can construct his or her own type of distance. The function
  method is more intuitive to work with than the older `makedist()`
  function.

- A new function,
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md),
  builds on the
  [`mdist()`](https://markmfredrickson.github.io/optmatch/dev/reference/mdist.md)
  structure to provide a convenient way to add calipers to a distance.
  In contrast to earlier ways of adding calipers,
  [`caliper()`](https://markmfredrickson.github.io/optmatch/dev/reference/caliper-methods.md)
  has an optional argument specify observations to be excluded from the
  caliper requirement — this permits one to relax it for just a few
  observations, for instance.

- [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  now removes strata in which matching failed (b/c the matching problem
  was found to be infeasible) before summarizing. It also indicates when
  such strata are present, and how many observations fall in them.

- Demo has been updated to reflect changes as of version 0.4, 0.5, 0.6.

### DEPRECATED & DEFUNCT

- The vignette is sufficiently out of date that it’s been removed.

### BUG FIXES

- subsetting of objects of class `Optmatch` now preserves
  matched.distances attribute.

- fixed bug in
  [`maxControlsCap()`](https://markmfredrickson.github.io/optmatch/dev/reference/minmaxctlcap.md)/[`minControlsCap()`](https://markmfredrickson.github.io/optmatch/dev/reference/minmaxctlcap.md)
  whereby they behaved unreliably on subclasses within which some
  subjects had no permissible matches.

- Removed unnecessary panic in
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  when it was given a `min.controls` argument with attributes other than
  names (as when it is created by
  [`tapply()`](https://rdrr.io/r/base/tapply.html)).

- fixed bug wherein
  [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  fails to retrieve balance tests if given a propensity model that had
  function calls in its formula.

- Documentation pages for
  [`fullmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
  [`pairmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md)
  filled out a bit.

## Changes in **optmatch** Version 0.5

### NEW FEATURES:

- [`summary.optmatch()`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)
  completely revised. It now reports information about the configuration
  of the matched sets and about matched distances. In addition, if given
  a fitted propensity model as a second argument it summarizes covariate
  balance.
