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
##* @param node_info NodeInfo specific to subproblem, or `NULL`
##* @return 
##* @keywords internal

solve_reg_fm_prob <- function(rownames, colnames, distspec, min.cpt,
                              max.cpt, tolerance, omit.fraction=NULL,
                              node_info = NULL)
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

  if (is.null(omit.fraction)) {
    f.ctls <- 1
  } else {
    if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction > 1) {
      stop("omit.fraction must be null or between 0 and 1")
    }
    f.ctls <- 1-omit.fraction
  }

    if (floor(min.cpt) > ceiling(max.cpt) |     #inconsistent max/min
        ceiling(1/min.cpt) < floor(1/max.cpt) | #controls per treatment
        !rfeas |  !cfeas  # either no controls or no treatments
        )
  {
    ans <- rep(NA_integer_,length(rownames)+length(colnames))
    names(ans) <- c(rownames, colnames)
    return(list(cells=ans, maxerr=NULL, distance=NULL))
  }
    

    old.o <- options(warn=-1)
    epsilon_lower_lim  <- max(dm$distance)/(.Machine$integer.max/64 -2)
    epsilon <- if (tolerance>0 & rfeas>1 & cfeas>1) {
                min(epsilon_lower_lim, tolerance/(rfeas + cfeas - 2))
            } else epsilon_lower_lim
    options(old.o)

    if (isTRUE(all.equal(dm[['distance']], 0)))
        dm[['distance']]  <- rep(1L, length(dm[['distance']])) # so we'll be routed to intSolve()

    temp <-
        if (is.integer(dm[['distance']]))
        {
            intSolve(dm, min.cpt, max.cpt, f.ctls, node_info)
        } else
        {
            doubleSolve(dm, min.cpt, max.cpt, f.ctls, node_info, rfeas, cfeas, epsilon)
        }

  temp$treated <- factor(temp$treated)
  temp$control <- factor(temp$control)
  ans <- rep(NA,length(rownames)+length(colnames))
  names(ans) <- c(rownames, colnames)

  matches <- solution2factor(temp)
  ans[names(matches)] <- matches

    if (!is.null(temp[["MCFSolution"]]))
        {
            temp[["MCFSolution"]]@subproblems[1L, "exceedance"]  <- temp$maxerr
            temp[["MCFSolution"]]@subproblems[1L, "feasible"]  <- any(temp$solutions==1L)

            ## Presently we can treat this subproblem as non-flipped even if it was,
            ## since `dm` will have been transposed in the event of flipping.  Doing
            ## so will prevent the evaluate_* routines below from getting confused.
            ## Of course we have to remember to set it based on actual information as
            ## the MCFsolutions object passes up through the point in the call stack
            ## where that transposition was made.
            temp[["MCFSolution"]]@subproblems[1L, "flipped"]  <- FALSE
            ## ... and now we can proceed with:
            evaluate_lagrangian(dm, temp[["MCFSolution"]]) ->
                temp[["MCFSolution"]]@subproblems[1L, "lagrangian_value"]
            evaluate_dual(dm, temp[["MCFSolution"]]) ->
                temp[["MCFSolution"]]@subproblems[1L,   "dual_value"    ]
            }

    return(list(cells = ans, err = temp$maxerr,
                MCFSolution=temp[["MCFSolution"]]
                )
           )
}


doubleSolve <- function(dm, min.cpt, max.cpt, f.ctls, node_info,
                        rfeas, cfeas, epsilon) 
{
    dm$distance  <- as.integer(ceiling(.5 + dm$distance / epsilon))
    if (!is.null(node_info))
        node_info$price  <- as.integer(ceiling(node_info$price / epsilon))
    
    intsol <- intSolve(dm=dm, min.cpt=min.cpt, max.cpt=max.cpt, f.ctls=f.ctls,
                       node_info = node_info)
    
    if (!is.null(intsol$MCFSolution))
    {
        intsol$MCFSolution@subproblems[1L,"resolution"]  <- epsilon
        
        intsol$MCFSolution@nodes[,'price']  <-
            intsol$MCFSolution@nodes[['price']] * epsilon
    }
    
    intsol$maxerr  <-
        if (any(is.na(intsol$solution))) { # i.e., problem was found infeasible.
            0 } else {
                  sum(intsol$solution * dm$distance, na.rm = TRUE) -
                      sum(intsol$solution * (intsol$distance - 1/2), na.rm = TRUE) * epsilon +
                      (sum(rfeas) > 1 & sum(cfeas) > 1) *
                      (sum(rfeas) + sum(cfeas) - 2 - sum(intsol$solution)) * epsilon
              }
    
  return(intsol)

}


intSolve <- function(dm, min.cpt, max.cpt, f.ctls, node_info = NULL)
    fmatch(dm, max.row.units = ceiling(1/min.cpt),
           max.col.units = ceiling(max.cpt),
           min.col.units = max(1, floor(min.cpt)),
           f=f.ctls, node_info =node_info)

##* Small helper function to turn a solution data.frame into a factor of matches
##* @keywords internal
solution2factor <- function(s) {
    s2  <- as.data.frame(s[c("control","treated","solution")])
  s2 <- s2[s2$solution == 1,]

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
