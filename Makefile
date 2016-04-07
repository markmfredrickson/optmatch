# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R
RCMD=R --vanilla -q -e

interactive:
	@$(LOAD) R -q --no-save

interactive-emacs:
	@$(LOAD) emacs -nw -f R

.devtools:
	@$(RCMD) "devtools:::$(FUNC)($(DEVTOOLSARG))"

DEVTOOLSARG=
dependencies: FUNC=install_deps
dependencies: DEVTOOLSARG=dependencies=TRUE
test: FUNC=test
check: FUNC=check
document: FUNC=document
vignette: FUNC=build_vignettes
clean-vignette: FUNC=clean_vignettes
build: FUNC=build
build_win: FUNC=build_win # Attempt a check & build on the win-builder server
dependencies test check document vignette clean-vignette build build_win: .devtools

clean: clean-vignette
	git clean -Xfd

spell-check-DESCRIPTION:
	aspell -c DESCRIPTION --personal=NULL
