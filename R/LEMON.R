#' (Internal) Helper function for accessing algorithms in LEMON solver
#'
#' @param algorithm LEMON algorithm to use. Choices are "CycleCancelling",
#'   "CapacityScaling", "CostScaling", "NetworkSimplex". Default is
#'   "NetworkSimplex".
#' @return String of the form "LEMON.<algorithm>"
#' @export
LEMON <- function(algorithm = "NetworkSimplex") {
  if (!(algorithm %in% c("CycleCancelling", "CapacityScaling", "CostScaling", "NetworkSimplex"))) {
    stop("Invalid LEMON algorithm. Valid algorithms are 'CycleCancelling', 'CapacityScaling', 'CostScaling', and 'NetworkSimplex'.")
  }
  return(paste0("LEMON.", algorithm))
}
