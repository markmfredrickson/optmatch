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
PKG=optmatch_$(VERSION)

# depend on the makefile so that updates to the version number will force a rebuild
# `git archive` doesn't export unarchived directories, so we export a .tar and untar it
$(PKG): Makefile
	rm -rf $(PKG)
	mkdir $(PKG)
	git archive --format=tar HEAD > $(PKG)/export.tar
	cd $(PKG) && tar xf export.tar
	rm $(PKG)/export.tar

$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION NAMESPACE ChangeLog NEWS R/* data/* demo/* inst/* man/* src/relax4s.f
	R --vanilla CMD Build $(PKG)

check: $(PKG).tar.gz
	R --vanilla CMD Check --no-multiarch $(PKG).tar.gz



