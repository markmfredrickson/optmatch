################################################################################
# Optmatch Class: the result of calling fullmatch()
################################################################################


####### Object Creation #########

makeOptmatch <- function(matching, # a list matches for different strata, each with names
                         treated, # list of the names of the treated units 
                         call) # the result of match.call() 
{    
  optmatch.obj <- unlist(mapply(function(label, groups) { paste(label, groups, sep = ".") }, 
                                   1:(length(matching)), 
                                   matching))
  optmatch.obj <- as.factor(optmatch.obj)

  names(optmatch.obj) <- names(unlist(matching))

  class(optmatch.obj) <- c("optmatch", "factor")

  # TODO: handle errors/failed matches
  # attr(strat.abv, "exceedances") <- err
  # if (sum(err, na.rm=TRUE)>TOL) {
  #     warning(
  #         paste("prescribed tol of ", tol, "per obs. poss. exceeded by up to ", 
  #           round(sum(err), 3), ".", sep="") )
  # }

  attr(optmatch.obj, "call") <- call
 
  attr(optmatch.obj, "contrast.group") <- names(optmatch.obj) %in% treated ### WHAT IS INROW?
  # TODO TURN ON WHEN MATCHED DISTANCES IS UPDATED
  # attr(optmatch.obj, "matched.distances") <- matched.distances(optmatch.obj, distance)
  
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

