#' @include feasible.R

.onAttach <- function(lib, pkg) {
packageStartupMessage("The optmatch package has an academic license. Enter relaxinfo() for more information.")
}

.onLoad <- function(lib, pkg) {
  setFeasibilityConstants()
  setTryRecovery()
  options("optmatch_verbose_messaging" = FALSE)
}
