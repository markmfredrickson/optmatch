pairmatch <- function(distance, controls = 1, remove.unmatchables = FALSE, ...) {

  # Should this checking be pushed to fullmatch to avoid duplication?
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification")
  }

  if (!all(floor(controls) == controls) | !all(controls > 0)) {
    stop("Minimum controls must be greater than treated units")  
  }

  subprobs <- findSubproblems(distance)

  if (length(controls) > 1 & !(length(subprobs) == length(controls))) {
    stop(paste("Controls argument must have same length as the number of subproblems (", 
      length(subprobs), ")", sep = ""))
  }

  omf <- mapply(controls, subprobs, FUN = function(control, prob) {
    # hard coding type based trimming for now. this should probably
    # be a DistanceSpecification method, e.g. finiteRows()
    if (remove.unmatchables) {
      if (inherits(prob, "matrix")) {
      # drop any rows that are entirely NA
      prob <- prob[apply(prob, 1, function(row) {
        any(is.finite(row)) }),]
      } else { 
        # assuming an InfinitySparseMatrix here      
        validrows <- which(1:(nrow(prob)) %in% prob@rows)    
        prob@dimension <- c(length(validrows), ncol(prob))
        prob@rownames <- prob@rownames[validrows]
      }
    }

    # a similar procedure is used to remove all control rows that
    # are unreachable

    if (inherits(prob, "matrix")) {
      # drop any rows that are entirely NA
      prob <- prob[, apply(prob, 2, function(row) {
        any(is.finite(row)) })]
    } else { 
        # assuming an InfinitySparseMatrix here      
        validcols <- which(1:(ncol(prob)) %in% prob@cols)    
        prob@dimension <- c(nrow(prob), length(validcols))
        prob@colnames <- prob@colnames[validcols]
    }

    nt <- nrow(prob)  
    nc <- ncol(prob)
    return((nc - control * nt)/nc)
  })

  if (any(omf<0)) {
    stop('not enough controls in some subclasses')
  }
  
  fullmatch(distance = distance, 
            min.controls = controls,
            max.controls = controls, 
            omit.fraction = omf,
            ...)
}

