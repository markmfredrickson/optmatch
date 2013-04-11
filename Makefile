################################################################################
# Development tools for optmatch
#
# The main targets available for use are: 
#
#   interactive: (the default) Launches an interactive R session with the current
#   working copy built into a package and automatically loaded.
#   
#   package: Builds a .tar.gz archive of the package suitable for installing,
#   etc.
#
#   spell: Checks the spelling of the package, see README.md for dependency
#   details.
#
#		test: Runs the 'testthat' package based tests and outputs any failing
#		tests. 'testthat' will be automatically installed in .local/ 
#
#		check: Performs `R CMD check` on the package, which is a slightly longer
#		process than just running the tests.
#
#		release: Builds the package, tests, and spellchecks in preparation for
#		sending to CRAN.
#
#		clean: removes all locally built files; to remove the downloaded
#		libraries, use `make clean-deps`.
#
# The version number of the package is set through the VERSION variable. This
# variable will change the name of the built .tar.gz file.
#
# All package dependencies, along with the optmatch package itself, will be
# installed in the .local directory and loaded throughout the build process.
# The default task (`make R`) also gives you access to these packages. If you
# want to test optmatch against a different version of a package, remove it
# from .local, install your version in .local, and finally touch the file
# ./local/<pkg>/INSTALLED. The last step is important to prevent your version
# of the package from being overwritten by the version from CRAN.
################################################################################
interactive: .local/optmatch/INSTALLED .local/testthat/INSTALLED .local/RItools/INSTALLED
	R_LIBS=.local R_PROFILE=load.R R -q --no-save 

### Package release scripts ###

VERSION=0.8-2
RELEASE_DATE=`date +%Y-%m-%d`
PKG=optmatch_$(VERSION)

# a useful helper for scripts who need to know what the package name is going to be
# use: R CMD INSTALL path/to/optmatch/$(cd path/to/optmatch && make current)
current: 
	@echo $(PKG).tar.gz

# depend on the makefile so that updates to the version number will force a rebuild
$(PKG): Makefile R/* tests/* inst/tests/* man/* inst/examples/*
	rm -rf $(PKG)
	rsync -a --exclude-from=.gitignore --exclude=.git* --exclude Makefile \
		--exclude=DESCRIPTION.template --exclude=NAMESPACE.static \
		--exclude=lexicon.txt --exclude=README.md --exclude=checkspelling.R \
		--exclude=optmatch.Rcheck \
		--exclude=vignettes \
		--exclude=load.R . $(PKG)

# You should probably use roxygen to add package dependecies, but if you must
# add them to DESCRIPTION.template
$(PKG)/DESCRIPTION: $(PKG) DESCRIPTION.template 
	sed s/VERSION/$(VERSION)/ DESCRIPTION.template | sed s/DATE/$(RELEASE_DATE)/ > $(PKG)/DESCRIPTION

# a macro for using the local directory only
LR = R_LIBS=.local R --vanilla

# Likewise, use roxygen to export functions
$(PKG)/NAMESPACE: $(PKG) $(PKG)/DESCRIPTION NAMESPACE.static .local/roxygen2/INSTALLED
	mkdir -p $(PKG)/man
	$(LR) -e "library(roxygen2); roxygenize('$(PKG)')"
	cat NAMESPACE.static >> $(PKG)/NAMESPACE

$(PKG).tar.gz: $(PKG) $(PKG)/DESCRIPTION $(PKG)/NAMESPACE NEWS R/* data/* demo/* inst/* man/* src/relax4s.f tests/*
	$(LR) CMD build $(PKG)

# a convenience target to get the current .tar.gz with having to know the
# version number in advance
package: $(PKG).tar.gz

# the spell task doesn't need the tar.gz particularly, but it does need DESCRIPTION and roxygen
spell: package 
	$(LR) -q --no-save -e "source('checkspelling.R') ; check_spelling('$(PKG)')"

lexicon.txt: package
	$(LR) -q --no-save -e "source('checkspelling.R') ; make_dictionary('$(PKG)')"

# the full (and slow) check process
check: $(PKG).tar.gz .local/testthat/INSTALLED .local/RItools/INSTALLED .local/biglm/INSTALLED
	$(LR) CMD check --as-cran --no-multiarch $(PKG).tar.gz

# getting ready to release
release: check spell
	git tag -a $(VERSION)
	@echo "Upload $(PKG).tar.gz to cran.r-project.org/incoming"
	@echo "Email to CRAN@R-project.org, subject: 'CRAN submission optmatch $(VERSION)'"

# depend on this file to decide if we need to install the local version
.local/optmatch/INSTALLED: $(PKG).tar.gz
	mkdir -p .local
	$(LR) CMD INSTALL --no-multiarch --library=.local $(PKG).tar.gz
	echo `date` > .local/optmatch/INSTALLED

# additional dependencies from CRAN
installpkg = mkdir -p .local ; $(LR) -e "install.packages('$(1)', repos = 'http://streaming.stat.iastate.edu/CRAN/')" ; date > .local/$(1)/INSTALLED
	
.local/devtools/INSTALLED:
	$(call installpkg,devtools)

.local/testthat/INSTALLED:
	$(call installpkg,testthat)

.local/RItools/INSTALLED:
	$(call installpkg,RItools)
	
.local/biglm/INSTALLED:
	$(call installpkg,biglm)

.local/profr/INSTALLED:
	$(call installpkg,profr)

# There is a bug in the released version of roxygen that prevents S4
# documentation from being properly built. This should be checked from time to
# time to see if the released version gets the bug fix.
.local/roxygen2/INSTALLED: .local/devtools/INSTALLED
	mkdir -p .local
	$(LR) -e "library(devtools) ; options(repos = 'http://streaming.stat.iastate.edu/CRAN/'); install_github(repo = 'roxygen', user = 'klutometis', branch = 's4',args=c('--no-multiarch'))"
	echo `date` > .local/roxygen2/INSTALLED

# test is just the internal tests, not the full R CMD Check
test: .local/optmatch/INSTALLED .local/testthat/INSTALLED .local/RItools/INSTALLED
	$(LR) -q -e "library(optmatch, lib.loc = '.local'); library(testthat); test_package('optmatch')"

# this will delete everything, except the CRAN dependencies in .local
clean:
	if [ -d .local ]; then mv .local .local-clean; fi
	git clean -Xfd
	if [ -d .local-clean ]; then mv .local-clean .local; fi

clean-deps:
	rm -rf .local

################################################################################
# Performance Testing
################################################################################

vignettes/performance/performance.pdf: .local/optmatch/INSTALLED .local/profr/INSTALLED \
																			 vignettes/performance/performance.Rnw \
																			 vignettes/performance/setup.rda \
																			 vignettes/performance/distance.rda \
																			 vignettes/performance/matching.rda \
																			 vignettes/performance/mdist.rda \
																			 vignettes/performance/scaling.rda 
	cd vignettes/performance && R_LIBS=../../.local R --vanilla CMD Sweave performance.Rnw
	cd vignettes/performance && latexmk -pdf performance.tex

vignettes/performance/setup.rda: .local/optmatch/INSTALLED vignettes/performance/setup.R
	cd vignettes/performance && R_LIBS=../../.local R --vanilla -f setup.R

vignettes/performance/distance.rda: vignettes/performance/setup.rda vignettes/performance/distance.R
	cd vignettes/performance && R_LIBS=../../.local R --vanilla -f distance.R

vignettes/performance/matching.rda: vignettes/performance/setup.rda vignettes/performance/matching.R
	cd vignettes/performance && R_LIBS=../../.local R --vanilla -f matching.R

vignettes/performance/mdist.rda: vignettes/performance/setup.rda vignettes/performance/mdist.R
	cd vignettes/performance && R_LIBS=../../.local R --vanilla -f mdist.R

vignettes/performance/scaling.rda: .local/optmatch/INSTALLED vignettes/performance/scaling.R
	cd vignettes/performance && R_LIBS=../../.local R --vanilla -f scaling.R
