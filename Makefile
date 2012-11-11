################################################################################
### Useful tasks for developing, not required to build the R package
################################################################################

R: .local/optmatch/INSTALLED
	R_PROFILE=load.R R -q --no-save 

### Package release scripts ###

VERSION=0.8-0
RELEASE_DATE=`date +%Y-%m-%d`
PKG=optmatch_$(VERSION)

# depend on the makefile so that updates to the version number will force a rebuild
# `git archive` doesn't export unarchived directories, so we export a .tar and untar it
# the code must be checked in to force a new export
$(PKG): Makefile R/* tests/* inst/tests/* man/* inst/examples/*
	rm -rf $(PKG)
	rsync -a --exclude-from=.gitignore --exclude=.git* --exclude Makefile \
		--exclude=DESCRIPTION.template --exclude=NAMESPACE.static \
		--exclude=lexicon.txt --exclude=README.md --exclude=checkspelling.R \
		--exclude=optmatch.Rcheck \
		--exclude=load.R . $(PKG)

$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

$(PKG)/NAMESPACE: $(PKG) $(PKG)/DESCRIPTION NAMESPACE.static .local/roxygen2/INSTALLED
	mkdir -p $(PKG)/man
	R_LIBS=.local R -e "library(roxygen2); roxygenize('$(PKG)')"
	cat NAMESPACE.static >> $(PKG)/NAMESPACE

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION $(PKG)/NAMESPACE NEWS R/* data/* demo/* inst/* man/* src/relax4s.f tests/*
	R --vanilla CMD build $(PKG)

package: $(PKG).tar.gz

# the spell task doesn't need the tar.gz particularly, but it does need DESCRIPTION and roxygen
spell: package 
	R -q --no-save -e "source('checkspelling.R') ; check_spelling('$(PKG)')"

lexicon.txt: package
	R -q --no-save -e "source('checkspelling.R') ; make_dictionary('$(PKG)')"

check: $(PKG).tar.gz 
	R --vanilla CMD Check --as-cran --no-multiarch $(PKG).tar.gz

release: check spell
	git tag -a $(VERSION)
	@echo "Upload $(PKG).tar.gz to cran.r-project.org/incoming"
	@echo "Email to CRAN@R-project.org, subject: 'CRAN submission optmatch $(VERSION)'"

# depend on this file to decide if we need to install the local version
.local/optmatch/INSTALLED: $(PKG).tar.gz
	mkdir -p .local
	R --vanilla CMD Install --no-multiarch --library=.local $(PKG).tar.gz
	echo `date` > .local/optmatch/INSTALLED

.local/roxygen2/INSTALLED:
	mkdir -p .local
	R_LIBS=.local R --vanilla -e "library(devtools) ; install_github(repo = 'roxygen', user = 'klutometis', branch = 's4',args=c('--no-multiarch'))"
	echo `date` > .local/roxygen2/INSTALLED

# test is just the internal tests, not the full R CMD Check
test: .local/optmatch/INSTALLED
	R --vanilla -q -e "library(optmatch, lib.loc = '.local'); library(testthat); test_package('optmatch')"

clean:
	git clean -xfd
