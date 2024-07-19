get_hints <- function(hint, subproblems)
{
  np <- length(subproblems)
  subproblemids <- names(subproblems)
  if (is.null(subproblemids)) subproblemids  <- character(1L)
  if (is.null(hint))
  {
    hints  <- rep(list(NULL), np)
  }
  else
  {
    hints  <- split(hint, hint[['groups']],
                    drop=TRUE # drops levels of hint$groups that aren't represented in hint
    )
    nohint  <- setdiff(subproblemids, names(hints))
    hints  <- hints[match(subproblemids, names(hints), 0L)]
    if (length(hints)>0) for (ii in 1L:length(hints)) hints[[ii]]  <- new("NodeInfo", hints[[ii]])
    if (length(nohint))
    {
      nullhint  <- rep(list(NULL), length(nohint))
      names(nullhint)  <- nohint
      hints  <- c(hints, nullhint)
      if (length(nohint)==np) warning("Hint lacks information about subproblems of this problem; ignoring.")
    }
    hints  <- hints[match(subproblemids, names(hints))]
  }
  hint.list <- mapply(prepare_subproblem_hint,
                      d = subproblems,
                      hint = hints,
                      SIMPLIFY = FALSE)
  return(hint.list)
}


get_tol_frac <- function(d, total_n, np)
{
  ncol <- dim(d)[2]
  nrow <- dim(d)[1]

  tol.frac <-
    if (total_n > 2 * np) {
      (nrow + ncol - 2)/(total_n - 2 * np)
    } else 1
  return(tol.frac)
}


get_flipped_status <- function(d, omf, mxctl, mnctl)
{
  ncol <- dim(d)[2]
  nrow <- dim(d)[1]

  if (switch(1 + is.na(omf), omf >= 0,  mxctl > .5)) {
    maxc <- min(mxctl, ncol)
    minc <- max(mnctl, 1/nrow)
    flipped  <- FALSE
  } else {
    maxc <- min(1/mnctl, ncol)
    minc <- max(1/mxctl, 1/nrow)
    d <- t(d)
    flipped  <- TRUE
  }
 return(list(subproblem = d,
             flipped_status = flipped,
             maxc = maxc,
             minc = minc))
}

prepare_subproblem_hint <- function(d, hint)
{
  if (is.null(hint)){
    hint  <- nodes_shell_fmatch(rownames(d),
                                colnames(d))
  }
  return(hint)
}

get_dms <- function(d,
                    hint)
{

  # if the subproblem is completely empty, short circuit
  if (length(d) == 0 ||
      all(is.infinite(d))) {
    return(NA)
  }

  rownames   <- subset(hint[["name"]], hint[['upstream_not_down']])
  colnames   <- subset(hint[["name"]], !hint[['upstream_not_down']])
  nrows  <- length(rownames)
  ncols  <- length(colnames)
  dm <- edgelist(d, c(rownames, colnames))
  return(dm)
}

get_epsilon <- function(subproblem, flipped,
                        hint, tolerance)
{
  rownames   <- subset(hint[["name"]], hint[['upstream_not_down']])
  colnames   <- subset(hint[["name"]], !hint[['upstream_not_down']])
  nrows  <- length(rownames)
  ncols  <- length(colnames)

  # matchable_nodes_info  <-
  #   filter(hint,
  #          is.na(hint[['upstream_not_down']]) | #retail bookkeeping nodes
  #            is_matchable(hint[['name']], dm, "either")
  #   )
  # rfeas  <- sum( matchable_nodes_info[['upstream_not_down']], na.rm=TRUE)
  # cfeas  <- sum(!matchable_nodes_info[['upstream_not_down']], na.rm=TRUE)
  #
  sub.dmat <- remove_inf_na_rows_cols(subproblem)
  if (!flipped)
  {
    rfeas <- nrow(sub.dmat)
    cfeas <- ncol(sub.dmat)
  } else {
    rfeas <- ncol(sub.dmat)
    cfeas <- nrow(sub.dmat)
  }


  max_dist <- max(sub.dmat, na.rm = TRUE)
  eps.val <- calculate_epsilon(rfeas,
                               cfeas,
                               tolerance,
                               max_dist)



  return(eps.val)
}

calculate_epsilon <- function(rfeas,
                              cfeas,
                              tolerance,
                              maxdist)
{
  old.o <- options(warn=-1)
  #epsilon_lower_lim  <- max(dm$'dist')/(.Machine$integer.max/64 -2)
  epsilon_lower_lim <- maxdist / (.Machine$integer.max/64 -2)
  epsilon <- if (tolerance>0 & rfeas>1 & cfeas>1) {
    max(epsilon_lower_lim, tolerance/(rfeas + cfeas - 2))
  } else epsilon_lower_lim
  options(old.o)
  return(epsilon)
}

parse_subproblems <- function(problems, min.controls,
                              max.controls, omit.fraction, hints,
                              total.n, TOL.in)
{
  np <- length(problems)

  #list of nodeInfo objects i think

  flipped_status <- mapply(get_flipped_status,
                           d = problems,
                           omf = omit.fraction,
                           mxctl = max.controls,
                           mnctl = min.controls,
                           SIMPLIFY = FALSE)

  problems <- lapply(flipped_status, function(x) x[['subproblem']])
  flippeds <- lapply(flipped_status, function(x) x[['flipped_status']])
  mincs <- lapply(flipped_status, function(x) x[['minc']])
  maxcs <- lapply(flipped_status, function(x) x[['maxc']])

  tol.fracs <- mapply(get_tol_frac,
                      d = problems,
                      total_n = rep(total.n, np),
                      np = rep(np, np))
  tolerances <- tol.fracs * TOL.in


  # subproblem.edgelists <- mapply(get_dms,
  #                     d = problems,
  #                     hint = hints,
  #                     SIMPLIFY = FALSE)
  # what about tolerances and epsilons when the problem is integer? account for this...

  epsilons <- mapply(get_epsilon,
                     hint = hints,
                     tolerance = tolerances,
                     subproblem = problems,
                     flipped = flippeds,
                     SIMPLIFY = FALSE)


  tmp <- list(subproblems = problems,
              flipped_status = flippeds,
              min.controls = mincs,
              max.controls = maxcs,
              epsilons = epsilons)
  return(tmp)
}


solve_reg_fm_prob2 <- function(node_info,
                              distspec,
                              min.cpt,
                              max.cpt,
                              omit.fraction=NULL,
                              solver,
                              epsilon) {

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


  if (all(abs(dm$'dist'- 0) < sqrt(.Machine$double.eps))) {
    dm$'dist'  <- rep(1L, length(dm$'dist')) # so we'll be routed to intSolve()
  }
  temp <-
    if (is.integer(dm$'dist')) {
      intSolve(dm, min.cpt, max.cpt, f.ctls, matchable_nodes_info, solver)
    } else {
      doubleSolve(dm, min.cpt, max.cpt, f.ctls, matchable_nodes_info,
                  rfeas, cfeas, epsilon, solver)
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



findSubproblemEpsilons <- function(tol = NULL,
                                   epsilon = NULL,
                                   subproblems,
                                   x,
                                   d)
{
  if (!is.null(tol)) #user specifies tolerance, rather than epsilon
  {
    total.n <- sum(dim(x))
    TOL <- tol * total.n
    np <- length(subproblems)
    ncol <- dim(d)[2]
    nrow <- dim(d)[1]


    tol.frac <-
      if (total.n > 2 * np) {
        (nrow + ncol - 2)/(total.n - 2 * np)
      } else 1

    # needs to handle multiple specifications for multiple subproblems, 1 specification for multiple subproblems, or 1 specification, 1 subproblem
    if (identical(length(tol), np)) {
      # if we have a tolerance specification for each subproblem
      tolerances <- TOL * tol.frac

    } else { #handle one tolerance for all subproblems. Need to check that other cases trigger errors
      if (length(tol) == 1)
      {

        tol.temp <- TOL * tol.frac #eventually, don't return tolerances, just return epsilons here
        tolerances <- rep(tol.temp, np)
      } else {
        stop("Please specify one tolerance for all subproblems, or one tolerance per subproblem.")
      }
    }
    return(tolerances)



  }

}

remove_inf_na_rows_cols <- function(mat) {
  # Identify rows that contain all Inf or NA
  rows_to_remove <- apply(mat, 1, function(row) all(is.infinite(row) | is.na(row)))

  # Identify columns that contain all Inf or NA
  cols_to_remove <- apply(mat, 2, function(col) all(is.infinite(col) | is.na(col)))

  # Remove identified rows and columns
  #cleaned_mat <- mat[!rows_to_remove, !cols_to_remove, drop = FALSE]
  cleaned_mat <- mat[!rows_to_remove, !cols_to_remove]
  return(cleaned_mat)
}
