##* Solves full matching problems that are regular,
##* in the sense that all row units are to be matched
##* but potentially some of the column units may be
##* left out.  Handles discretization of distances.
##* 
##* @param rownames character
##* @param colnames character
##* @param distspec InfinitySparseMatrix, matrix, etc (must have a `prepareMatching()` method)
##* @param min.cpt double, minimum permissible ratio of controls per treatment
##* @param max.cpt double, maximum permissible ratio of controls per treatment
##* @param tolerance 
##* @param omit.fraction 
##* @param matched.distances 
##* @param node.prices 
##* @param warm.start 
##* @param subproblemid 
##* @return 
##* @keywords internal

solve_reg_fm_prob <- function(rownames, colnames, distspec, min.cpt,
                         max.cpt, tolerance, omit.fraction=NULL, matched.distances=FALSE, node.prices = NULL, warm.start = NULL, subproblemid)
{

  if (min.cpt <=0 | max.cpt<=0) {
    stop("inputs min.cpt, max.cpt must be positive")
  }

  if (!all(rownames %in% dimnames(distspec)[[1]])) {
    stop("input \'rownames\' may only contain row names of input \'distspec\'")
  }

  if (!all(colnames %in% dimnames(distspec)[[2]])) {
    stop("input \'rownames\' may only contain col. names of input \'distspec\'")
  }

  # distance must have a prepareMatching object
  if (!hasMethod("prepareMatching", class(distspec))) {
    stop("Argument \'distspec\' must have a \'prepareMatching\' method")
  }

  # convert the distspec into a cannonical matching specification with columns
  # treated, control, distance

  dm <- prepareMatching(distspec)

  rownames <- as.character(rownames)
  colnames <- as.character(colnames)

  rfeas <- length(unique(dm$treated))
  cfeas <- length(unique(dm$control))
  # If any controls were unmatchable, they were dropped by prepareMatching, and
  # positive `omit.fraction`'s need to be updated.
  if (cfeas < length(colnames) & is.numeric(omit.fraction) && omit.fraction >0) {
    original_number_to_omit <- omit.fraction*length(colnames)
    number_implicitly_omitted_already <- length(colnames) - cfeas
    omit.fraction <- (original_number_to_omit - number_implicitly_omitted_already)/cfeas
    # This can happen if the number to be omitted is less than the number of unmatchables
    if (omit.fraction <= 0) {
      omit.fraction <- NULL
    }
  }

  # ... and similarly in the case of negative `omit.fraction` if there were
  # treatments that couldn't be matched.
  if (rfeas < length(rownames) & is.numeric(omit.fraction) && omit.fraction <0) {
    original_number_to_omit <- -1*omit.fraction*length(rownames)
    number_implicitly_omitted_already <- length(rownames) - rfeas
    omit.fraction <- - (original_number_to_omit - number_implicitly_omitted_already)/rfeas
    # This can happen if the number to be omitted is less than the number of unmatchables
    if (omit.fraction >= 0) {
      omit.fraction <- NULL
    }
  }

  if (floor(min.cpt) > ceiling(max.cpt) | ceiling(1/min.cpt) < floor(1/max.cpt))
  {
    ans <- rep("NA",length(rownames)+length(colnames))
    names(ans) <- c(rownames, colnames)
    return(list(cells=ans, maxerr=NULL, distance=NULL))
  }

  # the next block of code, the dm <- ... is commented out as
  # dm is no longer a matrix. Completely unreachable entries may be a
  # problem later, but
  if (is.null(omit.fraction)) {
    f.ctls <- 1
    # dm <- matrix(dm[rfeas, cfeas], sum(rfeas), sum(cfeas),
    # dimnames=list(rownames[rfeas], colnames[cfeas]))
  } else {
    if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction > 1) {
      stop("omit.fraction must be null or between 0 and 1")
    }

    f.ctls <- 1-omit.fraction
    # dm <- matrix(dm[rfeas,], sum(rfeas), length(colnames),
    #              dimnames=list(rownames[rfeas], colnames))
  }

  old.o <- options(warn=-1)
  if (any(dm$distance > 0)) {
    reso <- (.Machine$integer.max/64 -2)/max(dm$distance)
  } else {
    reso <- min(.Machine$integer.max/64 -2, (sum(rfeas)+sum(cfeas))/tolerance)
  }

  if (tolerance>0 & sum(rfeas)>1 & sum(cfeas)>1) {
    reso <- min(reso, (sum(rfeas) + sum(cfeas) - 2)/tolerance)
  }
    options(old.o)



  if (any(rfeas) & any(cfeas))
  {
    if(is.integer(dm[['distance']]) & any(dm$distance > 0)) #checking if all distances are integer
    {

      if(is.null(warm.start))
      {

        temp.with.nodes <- intSolve(dm, min.cpt, max.cpt, f.ctls, groupid = subproblemid)
      }
      else
      {
        temp.with.nodes <- intSolve(dm, min.cpt, max.cpt, f.ctls, node.prices, groupid = subproblemid)
      }

    }
    else
    {
      if(is.null(warm.start))
      {
        temp.with.nodes <- DoubleSolve(dm, rfeas, cfeas, min.cpt, max.cpt, tolerance, reso, f.ctls, groupid = subproblemid)
      }
      else
      {
        temp.with.nodes <- DoubleSolve(dm, rfeas, cfeas, min.cpt, max.cpt, tolerance, reso, f.ctls, warm.start = warm.start, groupid = subproblemid)
      }

    }



  } else {
    temp <- 0 ; maxerr <- 0 ; dist <- 0
  }
  temp <- temp.with.nodes$temp
  temp$treated <- factor(temp$treated)
  temp$control <- factor(temp$control)
  ans <- rep(NA,length(rownames)+length(colnames))
  names(ans) <- c(rownames, colnames)

  matches <- solution2factor(temp)
  ans[names(matches)] <- matches

    return(list(cells = ans, err = temp.with.nodes$maxerr#,
###                node.data = temp.with.nodes[["node.data"]],
###                prob.data = temp.with.nodes[["prob.data"]]
                )
           )
}


DoubleSolve <- function(dm, rfeas, cfeas, min.cpt,
                        max.cpt, tolerance, reso, f.ctls, warm.start = NULL, groupid = NULL) #warm.start should be a node.data data frame
{

  .matcher <- function(dm, toIntFunction, reso, min.cpt, max.cpt, f.ctls, nodeprices = NULL, groupid) {
    tmp <- dm
    tmp$distance <- toIntFunction(dm$distance * reso)

    if(!is.null(nodeprices))
    {
      node.ints <- toIntFunction(nodeprices * reso)
      obj <- intSolve(tmp, min.cpt, max.cpt, f.ctls, node.ints, groupid = groupid)
      return(obj)
    }
    else
    {
      obj <- intSolve(tmp, min.cpt, max.cpt, f.ctls, groupid = groupid)
      return(obj)
    }

  }

  # fmatch returns a matrix with columns `treatment`, `control`, and `solution`
  # it also has a column `distance` with toIntFuction(dm * reso)

  if(is.null(warm.start))
  {
    temp.with.nodes <- .matcher(dm, cadlag_ceiling, reso, min.cpt, max.cpt, f.ctls, groupid = groupid)
  }
  else
  {
    temp.with.nodes <- .matcher(dm, cadlag_ceiling, reso, min.cpt, max.cpt, f.ctls, nodeprices = warm.start, groupid = groupid)
  }


  if (any(is.na(temp.with.nodes$temp$solution))) {
    maxerr <- 0
  } else {
    maxerr <- sum(temp.with.nodes$temp$solution * dm$distance, na.rm = TRUE) -
      sum(temp.with.nodes$temp$solution * temp.with.nodes$temp$distance, na.rm = TRUE) / reso +
      (sum(rfeas) > 1 & sum(cfeas) > 1) *
      (sum(rfeas) + sum(cfeas) - 2 - sum(temp.with.nodes$temp$solution)) / reso
  }

##  if(!identical(options()$use_fallback_optmatch_solver, TRUE))
##  {
##    temp.with.nodes[["node.data"]]$price <- temp.with.nodes[["node.data"]]$price / reso
##  }

  temp.with.nodes$maxerr <- maxerr
##  temp.with.nodes[["prob.data"]]$tol = tolerance
##  temp.with.nodes[["prob.data"]]$reso = reso
##  temp.with.nodes[["prob.data"]]$exceedance <- maxerr
  return(temp.with.nodes)

}


intSolve <- function(dm, min.cpt, max.cpt, f.ctls, int.node.prices = NULL, groupid)
{
  if(!is.null(int.node.prices))
  {
    temp <- fmatch(dm, max.row.units = ceiling(1/min.cpt), max.col.units = ceiling(max.cpt), min.col.units = max(1, floor(min.cpt)), f=f.ctls, node_prices =int.node.prices)
  }
  else
  {
    temp <- fmatch(dm, max.row.units = ceiling(1/min.cpt), max.col.units = ceiling(max.cpt), min.col.units = max(1, floor(min.cpt)), f=f.ctls)
  }

###  temp.extended <- temp

  temp <- temp[1L:nrow(dm),] # Just the arcs representing potential matches

  match.with.node.prices <- list()
  match.with.node.prices[["temp"]] <- temp
###  if(identical(options()$use_fallback_optmatch_solver, TRUE))
###  {
###    match.with.node.prices[["node.data"]]
###  }
###  else
###  {
###    match.with.node.prices[["node.data"]] <-build_node_data(temp.extended = temp.extended, subproblemid = groupid)
###  }

  # not sure if following line should be one directly below this, or second option
###  match.with.node.prices[["prob.data"]] <- data.frame(max.control = max.cpt, min.control = min.cpt, omit.fraction = f.ctls, reso = NA, tol = NA, exceedance= 0, mean.control = NA, group = groupid)
###  #match.with.node.prices[["prob.data"]] <- data.frame(max.control = NA, min.control = NA, omit.fraction = NA, reso = NA, tol = NA, exceedance= 0, group = groupid)

  return(match.with.node.prices)
}


##' Rounding function like ceiling -- except at the integers themselves,
##' where it's always one ahead.  (Right instead of left continuous;
##' thus "cadlag".)  Point being to map nonnegative doubles
##' to _positive_ integers, as required by fmatch.
##'
##' The default tolerance value aims to steer well clear of inadvertently
##' passing 0L as where a 0 in x was being represented as a small negative
##' number, due to computer arithmetic. 
##' @title Right-continuous variant of base `ceiling()` function
##' @param x numeric
##' @param tol 
##' @return upwardly rounded version of x
##' @author Ben Hansen
##' @keywords internal
cadlag_ceiling  <- function(x, tol=.Machine$double.neg.eps^0.5) {1L + floor(x + tol)}


##* Small helper function to turn a solution data.frame into a factor of matches
##* @keywords internal
solution2factor <- function(s) {
  s2 <- s[s$solution == 1,]

  if (dim(s2)[1] == 0) {
    return(NULL)
  }

  # control units are labeled by the first treated unit to which they are connected
  # unlist(as.list(...)) was the best way I could find to make this into a vector, keeping names
  control.links <- unlist(as.list(by(s2, s2$control, function(x) { x[1,"treated"] })))

  # treated units are labeld by the label of the first control unit to which they are connected
  treated.links <- unlist(as.list(by(s2, s2$treated, function(x) { control.links[x[1, "control"]][1] })))

  # join the links
  return(c(treated.links, control.links))

}
