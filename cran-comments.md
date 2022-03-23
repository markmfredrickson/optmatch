# Test environments
* local OS X install, R 4.1.2
* win-builder (old, devel and release)
* mac-builer
* rhub

# R CMD check results

There were 2 note:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Josh Errickson <jerrick@umich.edu>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Klopfer (6:9)

CRAN repository db overrides:
  X-CRAN-Comment: Archived again on 2022-02-23 as still requires
    archived package 'RItools'.
  License_restricts_use: yes
  License_is_FOSS: no
Package license conflicts with 'License_is_FOSS: no' override
Package license conflicts with 'License_restricts_use: yes' override

Suggests or Enhances not in mainstream repositories:
  rrelaxiv
Availability using Additional_repositories specification:
  rrelaxiv   yes   https://josherrickson.github.io/rrelaxiv

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'rrelaxiv'
```

## Comments about NOTEs

This is a "new" release of an archived package.

"Klopfer" is a proper name.

See Submission Comments below; license is now standard and **rrelaxiv** is
non-CRAN package.

# Submission Comments

This is a re-release of the **optmatch** package which uses a standard license
as opposed to the non-standard license it was previously under.

This was accomplished by moving the code under the non-standard license into an
external package, **rrelaxiv**, which is listed under `Additional_repositories`
and will not be submitted to CRAN.

The package is fully functional without this additional package through the
addition of a dependence on the new **rlemon** package.
