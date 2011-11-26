matched.distances <- function(matchobj, distance, preserve.unit.names = FALSE)
{
  stopifnot(inherits(matchobj,"optmatch"))
  validDistanceSpecifcation(distance)

  # subsetting is a little rough right now, so making a temporary cast to matrix
  res <- tapply(names(matchobj), matchobj, FUN = function(x) {
      rs <- rownames(distance) %in% x
      cs <- colnames(distance) %in% x
      as.vector(subset(distance, subset = rs, select = cs, drop = !preserve.unit.names))
  })

  return(res)
}
