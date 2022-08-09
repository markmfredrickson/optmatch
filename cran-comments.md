# Test environments
* local OS X install, R 4.2.0
* win-builder (old, devel and release)
* mac-builer

# R CMD check results

There were 2 note:

```
❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Josh Errickson <jerrick@umich.edu>’

  Possibly misspelled words in DESCRIPTION:
    Klopfer (6:9)

❯ checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘RItools’
```

## Comments about NOTEs

Optional suggested dependency **RItools** is currently undergoing a revision to
return to CRAN. The main goal of this release is to ensure any use of RItools is
conditional so its absence of CRAN will not cause issues.

"Klopfer" is a proper name.

## Submission comments

This submission primarily addresses "Packages in Suggests should be used
conditionally".
