#' From a match (as produced by \code{pairmatch} or \code{fullmatch})
#' and a distance, extract the distances of matched units from their
#' matched counterparts.
#'
#' From a match (as produced by \code{pairmatch} or \code{fullmatch})
#' and a distance, extract the distances of matched units from their
#' matched counterparts.
#'
#' @title Determine distances between matched units
#' @param matchobj Value of a call to \code{pairmatch} or
#'   \code{fullmatch}.
#' @param distance Either a distance matrix or the value of a call to
#'   or \code{match_on}.
#' @param preserve.unit.names Logical.  If TRUE, for each matched set
#'   \code{matched.distances} returns the submatrix of the distance
#'   matrix corresponding to it; if FALSE, a vector containing the
#'   distances in that submatrix is returned.
#' @return A list of numeric vectors (or matrices) of distances, one
#'   for each matched set.  Note that a matched set with 1 treatment
#'   and k controls, or with k treatments and 1 control, has k, not
#'   k+1, distances.
#'
#' @author Ben B. Hansen
#' @examples
#'
#' data(plantdist)
#' plantsfm <- fullmatch(plantdist)
#' (plantsfm.d <- matched.distances(plantsfm,plantdist,pres=TRUE))
#' unlist(lapply(plantsfm.d, max))
#' mean(unlist(plantsfm.d))
#' @keywords nonparametric
#' @export
matched.distances <- function(matchobj, distance, preserve.unit.names = FALSE)
{
  stopifnot(inherits(matchobj,"optmatch"))
  validDistanceSpecification(distance) # this will stop() on invalid dist spec

  res <- tapply(names(matchobj), matchobj, FUN = function(x) {
      rs <- rownames(distance) %in% x
      cs <- colnames(distance) %in% x
      tmp <- as.vector(subset(distance, subset = rs, select = cs, drop = !preserve.unit.names))

      if (sum(rs) > sum(cs)) {
        names(tmp) <- rownames(distance)[rs]
      }

      if (sum(cs) > sum(rs)) {
        names(tmp) <- colnames(distance)[cs]
      }

      if (sum(cs) == sum(rs)) {
        names(tmp) <- colnames(distance)[cs] # either would be correct, and we don't make any guarantees about if we report t or c
      }
      return(tmp)
  })

  return(res)
}
