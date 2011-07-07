################################################################################
### Useful tasks for developing, not required to build the R package
################################################################################

local-install:
	rm -rf .local
	mkdir .local
	R CMD Install --library=.local .

autotest: local-install
	R -q -e "library(optmatch, lib.loc = '.local'); library(testthat); auto_test_package('.')"

build:
	R CMD Build .

check: build
	R CMD Check optmatch_0.7-2.tar.gz

clean:
	git clean

R: local-install
	R -q --no-save --interactive -e "library(optmatch, lib.loc = '.local')" 
