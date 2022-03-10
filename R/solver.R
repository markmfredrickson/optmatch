#' (Internal) Helper function for accessing algorithms in LEMON solver
#'
#' @param algorithm LEMON algorithm to use. Choices are "CycleCancelling",
#'   "CapacityScaling", "CostScaling", "NetworkSimplex". Default is
#'   "CycleCancelling".
#' @return String of the form "LEMON.<algorithm>"
#' @export
LEMON <- function(algorithm = "CycleCancelling") {
  if (!(algorithm %in% c("CycleCancelling", "CapacityScaling", "CostScaling", "NetworkSimplex"))) {
    stop("Invalid LEMON algorithm. Valid algorithms are 'CycleCancelling', 'CapacityScaling', 'CostScaling', and 'NetworkSimplex'.")
  }
  return(paste0("LEMON.", algorithm))
}

# Internal function to check the passed in solver.
#
# 1) If passed in "", check for presence of rrelaxiv. If present, return
# "RELAX-IV". If not present, return "LEMON.<defaultLEMON>" which is currently "LEMON.CycleCancelling".
# LEMON solver.
# 2) If passed in "RELAX-IV", check for presence of rrelaxiv. If present, return
# "RELAX-IV", else error.
# 3) If passed in "LEMON", return "LEMON.<defaultLEMON>" which is currently "LEMON.CycleCancelling".
# 4) If passed in LEMON(<solver>), evaluate and return "LEMON.<solver>" from LEMON().
# 5) All other input errors.
handleSolver <- function(solver) {
  hasrrelaxiv <- requireNamespace("rrelaxiv", quietly = TRUE)
  if (is.character(solver) & length(solver) == 1) {
    if (solver == "") {
      if (hasrrelaxiv) {
        return("RELAX-IV")
      } else {
        return(LEMON())
      }
    }

    if (solver == "RELAX-IV") {
      if (hasrrelaxiv) {
        return(solver)
      } else {
        stop("To use RELAX-IV solver, install package `rrelaxiv`.")
      }
    }

    if (solver == "LEMON") {
      return(LEMON())
    }

    if (solver %in% c("LEMON.CycleCancelling",
                      "LEMON.CostScaling",
                      "LEMON.CapacityScaling",
                      "LEMON.NetworkSimplex")) {
      return(solver)
    }

    stop("Invalid solver. Valid solvers are 'RELAX-IV' and 'LEMON'.")
  }
  stop("Invalid solver. Valid solvers are 'RELAX-IV' and 'LEMON'.")
}
