################################################################################
### Useful tasks for developing, not required to build the R package
################################################################################

autotest:
	rm -rf .autotest
	mkdir .autotest
	R CMD Install --library=.autotest .
	R -q -e "library(optmatch, lib.loc = '.autotest'); library(testthat); auto_test_package('.')"

