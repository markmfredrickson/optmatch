# Test environments
* local OS X install, R 4.4.1
* win-builder (old, devel and release)
* mac-builer

# R CMD check results

There were 2 note:

```

* checking CRAN incoming feasibility ... [16s] NOTE
Maintainer: 'Josh Errickson <jerrick@umich.edu>'

Suggests or Enhances not in mainstream repositories:
  rrelaxiv
Availability using Additional_repositories specification:
  rrelaxiv   yes   https://errickson.net/rrelaxiv
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'rrelaxiv'
```

The `rrelaxiv` package is available at the URL in DESCRIPTION, but is not on
CRAN as it contains a restrictive license.
