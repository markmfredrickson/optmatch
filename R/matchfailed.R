matchfailed <- function(matchobject) {
  grps <- attr(matchobject, "subproblem")
  failed <- sapply(split(matchobject, grps), function(x) { all(is.na(x)) })
  levels(grps) <- failed
  return(as.logical(grps))
} 
