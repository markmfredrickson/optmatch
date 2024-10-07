# Helper target for handling all devtools commands
.PHONY:.devtools
.devtools:
	@R -q -e "devtools:::$(FUNC)($(DEVTOOLSARG))"



# ******************************************************************************
# ******************** Interactive sessions ************************************
# ******************************************************************************

# For interacive sessions, load.R fixes a bug with devtool's `help` to enable
# `help` on functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R

.PHONY:interactive
interactive:
	@$(LOAD) R -q --no-save

.PHONY:interactive-emacs
interactive-emacs:
	@$(LOAD) emacs -nw -f R



# ******************************************************************************
# ******************** Working with Packages ***********************************
# ******************************************************************************

# Install dependencies
.PHONY:dependencies
dependencies: FUNC=install_deps
dependencies: DEVTOOLSARG=dependencies=TRUE

# Run all tests
.PHONY:test
test: FUNC=test

# Complete check
.PHONY:check
check: FUNC=check
check: DEVTOOLSARG=incoming=TRUE

# Complete check - no suggests included
.PHONY:check_no_suggests
check_no_suggests: FUNC=check
check_no_suggests: DEVTOOLSARG=incoming=TRUE,env_vars=c('_R_CHECK_DEPENDS_ONLY_'='true')

# Update documentation
.PHONY:document
document: FUNC=document

# Build vignettes
.PHONY:vignette
vignette: FUNC=build_vignettes

# Remove build vignettes
.PHONY:clean-vignette
clean-vignette: FUNC=clean_vignettes

# Build a source tarbell
.PHONY:build
build: FUNC=build
build: DEVTOOLSARG=args=c('--compact-vignettes=gs+qpdf')

# Update the pkgdown site
.PHONY:build_site
build_site: FUNC=build_site
build_site: DEVTOOLSARG=quiet=FALSE

dependencies test check check_no_suggests document: .devtools
vignette clean-vignette build build_site: .devtools



# ******************************************************************************
# ******************** Package Checks ******************************************
# ******************************************************************************

# Check on win-builder, old release
.PHONY:check_win_old
check_win_old: FUNC=check_win_oldrelease
check_win_old: DEVTOOLSARG=args=c('--compact-vignettes=gs+qpdf')

# Check on win-builder, current release
.PHONY:check_win
check_win: FUNC=check_win_release
check_win: DEVTOOLSARG=args=c('--compact-vignettes=gs+qpdf')

# Check on win-builder, development version
.PHONY:check_win_dev
check_win_dev: FUNC=check_win_devel
check_win_dev: DEVTOOLSARG=args=c('--compact-vignettes=gs+qpdf')

# Check on check_rhub
.PHONY:check_rhub
check_rhub: FUNC=check_rhub
check_rhub: DEVTOOLSARG=interactive=FALSE,env_vars=c('_R_CHECK_FORCE_SUGGESTS_'='false'),args=c('--compact-vignettes=gs+qpdf')

# Check on mac-builder (only release)
.PHONY:check_mac
check_mac: FUNC=check_mac_release
check_mac: DEVTOOLSARG=args=c('--compact-vignettes=gs+qpdf')

check_win check_win_dev check_win_old check_mac check_rhub: .devtools


# ******************************************************************************
# ******************** Release Package *****************************************
# ******************************************************************************

# Starts an interactive session which will take you through a collection of
# questions to ensure due-diliegnce, then submits to CRAN.

.PHONY:release
release:
	@echo Inside R, run "devtools::release(args=c('--compact-vignettes=gs+qpdf') )"
