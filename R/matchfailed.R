matchfailed <- function(match) {
  grps <- attr(match, "subproblem")
  failed <- sapply(split(match, grps), function(x) { all(is.na(x)) })
  levels(grps) <- failed
  return(as.logical(grps))
} 
