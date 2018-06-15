# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R
RCMD=R -q -e

.PHONY:interactive
interactive:
	@$(LOAD) R -q --no-save

.PHONY:interactive-emacs
interactive-emacs:
	@$(LOAD) emacs -nw -f R

.PHONY:.devtools
.devtools:
	@$(RCMD) "devtools:::$(FUNC)($(DEVTOOLSARG))"

DEVTOOLSARG=
.PHONY:dependencies
dependencies: FUNC=install_deps
dependencies: DEVTOOLSARG=dependencies=TRUE

.PHONY:test
test: FUNC=test

.PHONY:check
check: FUNC=check

.PHONY:document
document: FUNC=document

.PHONY:vignette
vignette: FUNC=build_vignettes

.PHONY:clean-vignette
clean-vignette: FUNC=clean_vignettes

.PHONY:build
build: FUNC=build

.PHONY:check_win_old
check_win_old: FUNC=check_win_oldrelease # Check & build on win-builder old release

.PHONY:check_win
check_win: FUNC=check_win_release        # ... on win-builder release

.PHONY:check_win_dev
check_win_dev: FUNC=check_win_devel    # ... on win-builder dev

dependencies test check document vignette clean-vignette build check_win check_win_dev check_win_oldrelease: .devtools

.PHONY:clean
clean: clean-vignette
	git clean -Xfd

.PHONY:spell-check-DESCRIPTION
spell-check-DESCRIPTION:
	aspell -c DESCRIPTION --personal=NULL
