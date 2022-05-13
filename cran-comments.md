# Test environments
* local OS X install, R 4.2.0
* win-builder (old, devel and release)
* mac-builer

# R CMD check results

There were 2 note:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Josh Errickson <jerrick@umich.edu>'

Suggests or Enhances not in mainstream repositories:
  rrelaxiv
Availability using Additional_repositories specification:
  rrelaxiv   yes   https://errickson.net/rrelaxiv

* checking package dependencies ... NOTE
Package suggested but not available for checking: 'rrelaxiv'
```

## Comments about NOTEs

Optional suggested dependency **rrelaxiv** has a very restrictive license and is
not being submitted to CRAN. It is listed under `Additional_repositories`.
