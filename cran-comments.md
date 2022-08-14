# Test environments
* local OS X install, R 4.2.1
* win-builder (old, devel and release)
* mac-builer

# R CMD check results

There were 2 note:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Josh Errickson <jerrick@umich.edu>'

Days since last update: 5

Suggests or Enhances not in mainstream repositories:
  rrelaxiv
Availability using Additional_repositories specification:
  rrelaxiv   yes   https://errickson.net/rrelaxiv
* checking package namespace information ... OK
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'rrelaxiv'
```

## Submission comments

I apologize for the quick re-release, but there's a failing test that was not
detected by any of my local checks, and was only found on a few CRAN checks.
This release fixes that failing test.
