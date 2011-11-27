matched.distances <- function(matchobj, distance, preserve.unit.names = FALSE)
{
  stopifnot(inherits(matchobj,"optmatch"))
  validDistanceSpecifcation(distance) # this will stop() on invalid dist spec

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
