################################################################################
# Optmatch Class: the result of calling fullmatch()
################################################################################


####### Object Creation #########

makeOptmatch <- function(distance, # a DistanceSpecification obj (e.g. a matrix or ISM)
                         solutions, # the results of calling SubDivStrat on multiple strata
                         call) # the result of match.call() 
{ 
  # pull out just the matching vectors
  matching <- lapply(solutions, function(x) { x$cells })

  treated <- rownames(distance)

  grpnames <- names(matching)
  if (is.null(grpnames)) {
    grpnames <- 1:(length(matching))  
  }
  
  optmatch.obj <- Reduce(mapply(function(label, groups) {
        tmp <- groups
        tmp[!is.na(groups)] <- paste(label, groups[!is.na(groups)], 
          sep = ".")
        return(tmp)
        }, grpnames, matching), f = c)

  optmatch.obj <- as.factor(optmatch.obj)
  names(optmatch.obj) <- unlist(sapply(matching, names)) 

  class(optmatch.obj) <- c("optmatch", "factor")

  tmp <- sapply(solutions, function(x) { x$err })
  names(tmp) <- grpnames
  attr(optmatch.obj, "exceedances") <- tmp

  attr(optmatch.obj, "call") <- call
 
  attr(optmatch.obj, "contrast.group") <- names(optmatch.obj) %in% treated ### WHAT IS INROW?
  # TODO TURN ON WHEN MATCHED DISTANCES IS UPDATED
  attr(optmatch.obj, "matched.distances") <- matched.distances(optmatch.obj, distance)
  
  return(optmatch.obj)
}


####### Subsetting and other manipulations #########

"[.optmatch" <-
  function(x, ..., drop=FALSE)
{
  y <- NextMethod("[")
  if  (!is.null(attr(x, "contrast.group"))) {
    cgs <- attr(x, "contrast.group")
    names(cgs) <- names(x)

    attr(y,"contrast.group") <- "["(cgs,...)
    names(attr(y, "contrast.group")) <-  NULL
  }

  ### The following is something of a kluge.  It would make more sense
  ### to remove matched distances that have been removed from the optmatch
  ### vector, but doing that is not straightforward, since the distances don't
  ### straighforwardly line up with the observations.  At present (version 0.6),
  ### the matched.distances attribute is only used in summary.optmatch;
  ### I have inserted code there to compensate for non-subsetting of the
  ### matched distances attribute in the case where matching has failed in some
  ### subclasses.
  if (!is.null(attr(x, "matched.distances"))) {
    attr(y, "matched.distances") <- attr(x, "matched.distances")
  }
  
  class(y) <- c("optmatch", "factor")
  
  return(y)
}

