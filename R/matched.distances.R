matched.distances <- function(matchobj, distance, preserve.unit.names = FALSE)
{
  stopifnot(inherits(matchobj,"optmatch"))
  validDistanceSpecifcation(distance) # this will stop() on invalid dist spec

  res <- tapply(names(matchobj), matchobj, FUN = function(x) {
      rs <- rownames(distance) %in% x
      cs <- colnames(distance) %in% x
      tmp <- as.vector(subset(distance, subset = rs, select = cs, drop = !preserve.unit.names))
      return(tmp)
  })

  return(res)
}
