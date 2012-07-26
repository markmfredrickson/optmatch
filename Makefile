################################################################################
### Useful tasks for developing, automating release process
################################################################################

### Interactive use and testing targets ###

R: OPTMATCH_TIMESTAMP
	R_PROFILE=interactive.R R -q --no-save 
	
.local:
	mkdir .local

OPTMATCH_TIMESTAMP: .local R/* tests/* inst/tests/*
	R --vanilla CMD Install --no-multiarch --library=.local .
	date > OPTMATCH_TIMESTAMP

test: OPTMATCH_TIMESTAMP
	R --vanilla -q -e "library(optmatch, lib.loc = '.local'); library(testthat); test_package('optmatch')"

clean:
	git clean

### Package release scripts ###

VERSION=0.7-3
RELEASE_DATE=`date +%Y-%m-%d`

# depend on the makefile so that updates to the version number will force a rebuild of DESCRIPTION
DESCRIPTION: DESCRIPTION.template Makefile
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > DESCRIPTION

optmatch_$(VERSION).tar.gz: DESCRIPTION NAMESPACE ChangeLog NEWS R/* data/* demo/* inst/* man/* src/relax4s.f
	R --vanilla CMD Build .

check: optmatch_$(VERSION).tar.gz
	R --vanilla CMD Check --no-multiarch optmatch_$(VERSION).tar.gz



