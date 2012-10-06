matchfailed <- function(x) {
  grps <- attr(x, "subproblem")
  failed <- sapply(split(x, grps), function(x) { all(is.na(x)) })
  levels(grps) <- failed
  return(as.logical(grps))
} 
