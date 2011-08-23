################################################################################
### Useful tasks for developing, not required to build the R package
################################################################################

R: OPTMATCH_TIMESTAMP
	R -q --no-save 
	
.local:
	mkdir .local

OPTMATCH_TIMESTAMP: .local R/* tests/* inst/tests/*
	R --vanilla CMD Install --library=.local .
	date > OPTMATCH_TIMESTAMP

autotest: OPTMATCH_TIMESTAMP
	R --vanilla -q -e "library(optmatch, lib.loc = '.local'); library(testthat); auto_test('./R', './inst/tests', 'summary')"

build:
	R --vanilla CMD Build .

check: build
	R --vanilla CMD Check optmatch_0.7-2.tar.gz

clean:
	git clean


