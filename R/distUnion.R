#' Combine multiple distance specifications into a single distance specification.
#'
#' For combining multiple distance specifications with common
#' controls, but different treated units, \code{\link{rbind}} provides
#' a way to combine the different objects. Likewise,
#' \code{\link{cbind}} provides a way to combine distance
#' specifications over common treated units, but different control
#' units.
#'
#' \code{distUnion} can combine distance units that have common
#' treated and control units into a coherent single distance object.
#' 
#' @param ... The distance specifications.
#' @return An InfinitySparseMatrix object with all treated and control
#' units from the arguments.
#' @export
distUnion <- function(...) {
  arglst <- list(...)

  if (!all(sapply(arglst, validDistanceSpecification))) {
    stop("All arguments must be valid distance specifications")
  }

  isms <- lapply(arglst, as.InfinitySparseMatrix) 

  treateds <- lapply(isms, function(i) { i@rownames })
  controls <- lapply(isms, function(i) { i@colnames })

  utreated  <- unique(unlist(treateds))
  ucontrols <- unique(unlist(controls))

  tmap <- 1:length(utreated)
  names(tmap) <- utreated

  cmap <- 1:length(ucontrols)
  names(cmap) <- ucontrols

  updated.isms <- lapply(isms, function(i) {
    rnms <- i@rownames[i@rows]
    cnms <- i@colnames[i@cols]
    i@rows <- tmap[rnms]
    i@cols <- cmap[cnms]
    return(i)
  })

  unionism <- makeInfinitySparseMatrix(
      unlist(sapply(updated.isms, function(i) { i@.Data })),
      unlist(sapply(updated.isms, function(i) { i@cols })),
      unlist(sapply(updated.isms, function(i) { i@rows })),
      ucontrols,
      utreated)
  
  return(unionism)
}
