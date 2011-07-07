pairmatch <- function(distance, controls=1, tol=0.001, remove.unmatchables=FALSE) {

  # Should this checking be pushed to fullmatch to avoid duplication?
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification")
  }

  if (!all(floor(controls) == controls) | !all(controls > 0)) {
    stop("Minimum controls must be greater than treated units")  
  }

n.Tx <- if (remove.unmatchables)
  function(x) {sum(apply(x,1,function(row) any(is.finite(row))))} else nrow 
                                 
if (inherits(distance, "optmatch.dlist"))
  {
nt <- sapply(distance, n.Tx)
nc <- sapply(distance, function(x){ncol(x)})
omf <- (nc-controls*nt)/nc
if (any(omf<0)) stop('not enough controls in some subclasses')
} else {
  nt <- n.Tx(distance)
  nc <- ncol(distance)
  omf <- (nc-controls*nt)/nc
  if (any(omf<0)) stop('not enough controls')
}
fullmatch(distance=distance, min.controls=controls,
          max.controls=controls, omit.fraction=omf,
          tol=tol)
  }

