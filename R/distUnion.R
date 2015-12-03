#' Combine multiple distance specifications into a single distance specification.
#'
#' Creates a new distance specification from the union of two or more
#' distance specifications. The constituent distances specifications
#' may have overlapping treated and control units (identified by the
#' \code{rownames} and \code{colnames} respectively).
#'
#' For combining multiple distance specifications with common
#' controls, but different treated units, \code{\link{rbind}} provides
#' a way to combine the different objects. Likewise,
#' \code{\link{cbind}} provides a way to combine distance
#' specifications over common treated units, but different control
#' units.
#'
#' \code{distUnion} can combine distance units that have common
#' treated and control units into a coherent single distance
#' object. If there are duplicate treated-control entries in multiple
#' input distances, the first entry will be used.
#' 
#' @param ... The distance specifications (as created with with
#' \code{\link{match_on}}, \code{\link{exactMatch}}, or other distance
#' creation function).
#' @return An InfinitySparseMatrix object with all treated and control
#' units from the arguments combined. Duplicate entries are resolved
#' in favor of the earliest argument (e.g., \code{distUnion(A, B)}
#' will favor entries in \code{A} over entries in \code{B}).
#' @seealso \code{\link{match_on}}, \code{\link{exactMatch}},
#' \code{\link{fullmatch}}, \code{\link{pairmatch}},
#' \code{\link{cbind}}, \code{\link{rbind}}
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


  pairs <- matrix(c(unlist(sapply(updated.isms, function(i) { i@cols })),
                    unlist(sapply(updated.isms, function(i) { i@rows }))),
                  ncol = 2)

  dups <- duplicated(pairs)

  pairs <- pairs[!dups, ]

  unionism <- makeInfinitySparseMatrix(
      unlist(sapply(updated.isms, function(i) { i@.Data }))[!dups],
      cols = pairs[, 1],
      rows = pairs[, 2],
      ucontrols,
      utreated)
  
  return(unionism)
}
