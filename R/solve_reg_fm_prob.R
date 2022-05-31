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
##* @param node_info NodeInfo. (Only names and upstream_not_down cols will be used.)
##* @param distspec InfinitySparseMatrix, matrix, etc (must have a `prepareMatching()` method)
##* @param min.cpt double, minimum permissible ratio of controls per treatment
##* @param max.cpt double, maximum permissible ratio of controls per treatment
##* @param tolerance 
##* @param omit.fraction 
##* @return 
##* @keywords internal

solve_reg_fm_prob <- function(node_info, distspec, min.cpt,
                              max.cpt, tolerance, omit.fraction=NULL
                              )
{

  if (min.cpt <=0 | max.cpt<=0) {
    stop("inputs min.cpt, max.cpt must be positive")
  }
    stopifnot(is(node_info, "NodeInfo"))
    rownames   <- subset(node_info[["name"]], node_info[['upstream_not_down']])
    colnames   <- subset(node_info[["name"]], !node_info[['upstream_not_down']])
    nrows  <- length(rownames)
    ncols  <- length(colnames)
  if (!all(rownames %in% dimnames(distspec)[[1]])) {
    stop("node_info rownames must be rownames for \'distspec\'")
  }

  if (!all(colnames %in% dimnames(distspec)[[2]])) {
    stop("node_info colnames must be colnames for \'distspec\'")
  }
 
  # distance must have an edgelist method
  if (!hasMethod("edgelist", class(distspec))) {
    stop("Argument \'distspec\' must have a \'edgelist\' method")
  }

  # convert the distspec to cannonical EdgeList
  dm <- edgelist(distspec, c(rownames, colnames))


  matchable_nodes_info  <-
      filter(node_info,
             is.na(node_info[['upstream_not_down']]) | #retail bookkeeping nodes
             is_matchable(node_info[['name']], dm, "either")
             )
  rfeas  <- sum( matchable_nodes_info[['upstream_not_down']], na.rm=TRUE)
  cfeas  <- sum(!matchable_nodes_info[['upstream_not_down']], na.rm=TRUE)
  ## Update `omit.fraction`:
  if (cfeas < ncols & is.numeric(omit.fraction) && omit.fraction >0) {
      original_number_to_omit <- omit.fraction*ncols
      number_implicitly_omitted_already <- ncols - cfeas
      omit.fraction <-
          (original_number_to_omit - number_implicitly_omitted_already)/cfeas
      ## If the number to be omitted is less than the number of unmatchable
      ## columns, ncols - cfeas, then this omit.fraction can be negative.  Then
      ## we just omit a few more than directed by the supplied omit.fraction:
      if (omit.fraction <= 0) omit.fraction <- NULL
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
    ans <- rep(NA_integer_, nrows + ncols)
    names(ans) <- c(rownames, colnames)
    return(list(cells=ans, err=0, MCFSolution=NULL))
  }
    

    old.o <- options(warn=-1)
    epsilon_lower_lim  <- max(dm$'dist')/(.Machine$integer.max/64 -2)
    epsilon <- if (tolerance>0 & rfeas>1 & cfeas>1) {
                min(epsilon_lower_lim, tolerance/(rfeas + cfeas - 2))
            } else epsilon_lower_lim
    options(old.o)

    if (all(abs(dm$'dist'- 0) < sqrt(.Machine$double.eps)))
    {
        dm$'dist'  <- rep(1L, length(dm$'dist')) # so we'll be routed to intSolve()
        }
    temp <-
        if (is.integer(dm$'dist'))
        {
            intSolve(dm, min.cpt, max.cpt, f.ctls, matchable_nodes_info)
        } else
        {
            doubleSolve(dm, min.cpt, max.cpt, f.ctls, matchable_nodes_info, rfeas, cfeas, epsilon)
        }

  temp$treated <- factor(temp[['i']]) # levels of these factors now convey
  temp$control <- factor(temp[['j']]) # treatment/control distinction
  matches <- solution2factor(temp)

  ans <- rep(NA_character_, nrows + ncols)
  names(ans) <- c(rownames, colnames)
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
            nodeinfo(temp[["MCFSolution"]])  <-
                update(node_info, nodeinfo(temp[["MCFSolution"]]))
            }

    return(list(cells = ans, err = temp$maxerr,
                MCFSolution=temp[["MCFSolution"]]
                )
           )
}


doubleSolve <- function(dm, min.cpt, max.cpt, f.ctls, node_info,
                        rfeas, cfeas, epsilon) 
{
    dm_distance  <- dm$'dist' #Used below in roundoff error crude estimate
    dm$'dist'  <- as.integer(ceiling(.5 + dm$'dist' / epsilon))
    if (!is.null(node_info))
        node_info$price  <-
            as.integer(ifelse(abs(node_info$price) <sqrt(.Machine$double.eps),
                              0, 
                              ceiling(node_info$price / epsilon)
                              )
                       )
    
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
            0 } else { #roundoff error crude estimate
                  sum(intsol$solution * dm_distance, na.rm = TRUE) -
                      sum(intsol$solution * (intsol$dist - 1/2), na.rm = TRUE) * epsilon +
                      (sum(rfeas) > 1 & sum(cfeas) > 1) *
                      (sum(rfeas) + sum(cfeas) - 2 - sum(intsol$solution)) * epsilon
              }
    
  return(intsol)

}


intSolve <- function(dm, min.cpt, max.cpt, f.ctls, node_info)
{
    stopifnot(is(dm, "EdgeList"), is(node_info, "NodeInfo"))
    if (!is.integer(node_info[['price']]))
    {
        price_col_position  <- which(node_info@names=="price")
        node_info@.Data[[price_col_position]]  <-
            as.integer(round(node_info[['price']]))
    }

    fmatch(dm, max.row.units = ceiling(1/min.cpt),
           max.col.units = ceiling(max.cpt),
           min.col.units = max(1, floor(min.cpt)),
           f=f.ctls, node_info =node_info)
    }

##* Small helper function to turn a solution data.frame into a factor of matches
##* @keywords internal
solution2factor <- function(s) {
    s2  <- as.data.frame(s[c("control","treated","solution")])
  s2 <- s2[s2[['solution']] == 1,,drop=FALSE]

  if (dim(s2)[1] == 0) {
    return(NULL)
  }

  ## control units are labeled by the first treated unit
  ## to which they are connected.
  control.links  <- factor(character(nlevels(s2[['control']])),
                           levels=levels(s2[['treated']])
                           )
  names(control.links)  <- levels(s2[['control']])
  reduced  <- s2[!duplicated(s2[['control']]),,drop=FALSE]
  control.links[as.character(reduced[['control']])]  <- 
      reduced[['treated']]

  ## Treated units are labeled the same as the
  ## controls they're connected to.  When
  ## a treatment is matched to multiple controls,
  ## they each have the same label, so we just
  ## use the first instance.
  treated.links  <- factor(character(nlevels(s2[['treated']])),
                           levels=levels(s2[['treated']])
                           )
  names(treated.links)  <- levels(s2[['treated']])
  reduced  <- s2[!duplicated(s2[['treated']]), , drop=FALSE]
  treated.links[as.character(reduced[['treated']])]  <- 
    control.links[as.character(reduced[['control']])]

  ## join the links. (Note that `c()` drops the factor
  ## structure, giving a named integer vector.)
  c(treated.links, control.links)
}

##* Make shell of node table, as required by `fmatch()`
##*
##* For now, only name and upstream_not_down cols are meaningful
##* @title MCF node table for ordinary full matches 
##* @param rownames character
##* @param colnames character
##* @return NodeInfo
##* @author Hansen
nodes_shell_fmatch <- function(rownames, colnames) {
    dm  <- c(length(rownames), length(colnames))
    ans  <- data.frame(name=c(rownames, colnames,"(_Sink_)", "(_End_)"),
                       price=0,
                       upstream_not_down=c(rep(c(TRUE, FALSE),
                                               times=dm),
                                           rep(NA, 2)
                                           ),
                       supply=integer(sum(dm)+2),
                       stringsAsFactors=FALSE
                       )
    new("NodeInfo", ans)
          }
          

