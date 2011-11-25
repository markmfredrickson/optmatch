matched.distances <- function(matchobj, distance, preserve.unit.names = FALSE)
{
  stopifnot(inherits(matchobj,"optmatch"))
  validDistanceSpecifcation(distance)

  # subsetting is a little rough right now, so making a temporary cast to matrix
  distance <- as.matrix(distance)
  res <- tapply(names(matchobj), matchobj, FUN = function(x) {
      distance[match(x, dimnames(distance)[[1]], nomatch = 0),
               match(x, dimnames(distance)[[2]], nomatch = 0),
               drop = !preserve.unit.names]
  })

  return(res)
}
