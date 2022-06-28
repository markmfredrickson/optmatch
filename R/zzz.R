#' @include feasible.R

.onLoad <- function(lib, pkg) {
  setFeasibilityConstants()
  setTryRecovery()
  options("optmatch_verbose_messaging" = FALSE)
}
