pairmatch <- function(distance, controls=1, tol=0.001, remove.unmatchables=FALSE) {

  # Should this checking be pushed to fullmatch to avoid duplication?
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification")
  }

  if (!all(floor(controls) == controls) | !all(controls > 0)) {
    stop("Minimum controls must be greater than treated units")  
  }

  # hard coding type based trimming for now. this should probably
  # be a DistanceSpecification method, e.g. finiteRows()
  if (remove.unmatchables) {
    if (inherits(distance, "matrix")) {
      # drop any rows that are entirely NA
      distance <- distance[apply(distance, 1, function(row) {
        any(is.finite(row)) }),]
    } else { 
      # assuming an InfinitySparseMatrix here      
      validrows <- which(1:(nrow(distance)) %in% distance@rows)    
      distance@dimension <- c(length(validrows), ncol(distance))
      distance@rownames <- distance@rownames[validrows]
    }
  }

  nt <- nrow(distance)  
  nc <- ncol(distance)
  omf <- (nc - controls * nt)/nc
  if (any(omf<0)) {
    stop('not enough controls in some subclasses')
  }
  
  fullmatch(distance=distance, min.controls=controls,
    max.controls=controls, omit.fraction=omf,
    tol=tol)
}

