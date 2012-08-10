################################################################################
### Useful tasks for developing, automating release process
################################################################################

R: .local/optmatch/INSTALLED
	R_PROFILE=interactive.R R -q --no-save 

### Package release scripts ###

VERSION=0.7-6
RELEASE_DATE=`date +%Y-%m-%d`
PKG=optmatch_$(VERSION)

# depend on the makefile so that updates to the version number will force a rebuild
# `git archive` doesn't export unarchived directories, so we export a .tar and untar it
# the code must be checked in to force a new export
$(PKG): Makefile R/* tests/*  man/*
	rm -rf $(PKG)
	rsync -a --exclude-from=.gitignore --exclude=.git* --exclude Makefile --exclude=DESCRIPTION.template --exclude=interactive.R . $(PKG)

$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION NAMESPACE ChangeLog NEWS R/* data/* demo/* inst/* man/* src/relax4s.f tests/*
	R --vanilla CMD Build $(PKG)

check: $(PKG).tar.gz
	R --vanilla CMD Check --as-cran --no-multiarch $(PKG).tar.gz

release: check
	git tag -a $(VERSION)
	@echo "Upload $(PKG) to cran.r-project.org/incoming"
	@echo "Email to CRAN@R-project.org, subject: 'CRAN submission optmatch $(VERSION)'"

# depend on this file to decide if we need to install the local version
.local/optmatch/INSTALLED: $(PKG).tar.gz
	mkdir -p .local
	R --vanilla CMD Install --no-multiarch --library=.local $(PKG).tar.gz
	echo `date` > .local/optmatch/INSTALLED

# test is just the internal tests, not the full R CMD Check
test: .local/optmatch/INSTALLED
	R --vanilla -q -e "library(optmatch, lib.loc = '.local'); library(testthat); test_package('optmatch')"

clean:
	git clean -xfd
