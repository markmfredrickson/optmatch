##* Solves full matching problems that are regular,
##* in the sense that all row units are to be matched
##* but potentially some of the column units may be
##* left out.
##*
##* Structured so that doubleSolve calls IntSolve after
##* problem has been converted to integer resolution, and
##* then IntSolve calls fmatch.
##* This function handles:
##* \itemize{
##* \item discretization of distances, as needed, along
##*       with conversion of allocated regret budget into
##8       a resolution for the discretization;
##* \item reconciliation of omit.fraction w/ necessary
##*       exclusions due to unit-wise unmatchability;
##* \item flagging of problems with constraints that
##*       can be seen to be infeasible even w/o calling
##*       the solver.
##* }
##* @param rownames character
##* @param colnames character
##* @param distspec InfinitySparseMatrix, matrix, etc (must have a `prepareMatching()` method)
##* @param min.cpt double, minimum permissible ratio of controls per treatment
##* @param max.cpt double, maximum permissible ratio of controls per treatment
##* @param tolerance 
##* @param omit.fraction 
##* @param matched.distances 
##* @param warm.start Numeric vector of node prices
##* @param subproblemid 
##* @return 
##* @keywords internal

solve_reg_fm_prob <- function(rownames, colnames, distspec, min.cpt,
                              max.cpt, tolerance, omit.fraction=NULL, matched.distances=FALSE,
                              warm.start = NULL, subproblemid)
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

    if (floor(min.cpt) > ceiling(max.cpt) | ceiling(1/min.cpt) < floor(1/max.cpt) |
        !rfeas |  !cfeas )
  {
    ans <- rep("NA",length(rownames)+length(colnames))
    names(ans) <- c(rownames, colnames)
    return(list(cells=ans, maxerr=NULL, distance=NULL))
  }

  if (is.null(omit.fraction)) {
    f.ctls <- 1
  } else {
    if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction > 1) {
      stop("omit.fraction must be null or between 0 and 1")
    }
    f.ctls <- 1-omit.fraction
  }

  old.o <- options(warn=-1)
  if (any(dm$distance > 0)) {
    reso <- (.Machine$integer.max/64 -2)/max(dm$distance)
  } else {
    reso <- min(.Machine$integer.max/64 -2, (rfeas+cfeas)/tolerance)
  }

  if (tolerance>0 & rfeas>1 & cfeas>1) {
    reso <- min(reso, (rfeas + cfeas - 2)/tolerance)
  }
  options(old.o)



    temp.with.nodes <-
        if (is.integer(dm[['distance']]) & any(dm$distance > 0)) #checking if all distances are integer
        {
            intSolve(dm, min.cpt, max.cpt, f.ctls, int.node.prices=warm.start, groupid = subproblemid)
        } else
        {
            doubleSolve(dm, rfeas, cfeas, min.cpt, max.cpt, tolerance, reso, f.ctls, warm.start = warm.start, groupid = subproblemid)
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


doubleSolve <- function(dm, rfeas, cfeas, min.cpt,
                        max.cpt, tolerance, reso, f.ctls, warm.start = NULL, groupid = NULL) #warm.start should be a node.data data frame
{

    dm$distance  <- cadlag_ceiling(dm$distance * reso)
    node.ints  <- if (is.null(warm.start)) NULL else cadlag_ceiling(warm.start * reso)
    
    temp.with.nodes <- intSolve(dm=dm, min.cpt=min.cpt, max.cpt=max.cpt, f.ctls=f.ctls,
                                int.node.prices = node.ints, groupid = groupid)


  if (any(is.na(temp.with.nodes$temp$solution))) { # i.e., problem was found infeasible.
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
    temp <- fmatch(dm, max.row.units = ceiling(1/min.cpt), max.col.units = ceiling(max.cpt), min.col.units = max(1, floor(min.cpt)), f=f.ctls, node_prices =int.node.prices)

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
