#' @include feasible.R

.onAttach <- function(lib, pkg) {
packageStartupMessage(paste(
"You're loading optmatch, by B. Hansen and M. Fredrickson.\n",
"The optmatch package makes essential use of D. P. Bertsekas\n",
"and P. Tseng\'s RELAX-IV algorithm and code, as well as\n" ,
"Bertsekas\' AUCTION algorithm and code.  Using the software\n",
"to \'satisfy in any part commercial delivery requirements to\n",
"government or industry\' requires a special agreement with\n",
"Dr. Bertsekas. For more information, enter\n",
"relaxinfo() at the command line.\n"
))
}

.onLoad <- function(lib, pkg) {
  setFeasibilityConstants()
  setTryRecovery()
  options("optmatch_verbose_messaging" = FALSE)
}
