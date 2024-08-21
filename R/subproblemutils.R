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

  sub.dmat <- as.InfinitySparseMatrix(subproblem)

  if (!flipped)
  {
    rfeas <- length(unique(sub.dmat@rows))
    cfeas <- length(unique(sub.dmat@cols))
  } else {
    rfeas <- length(unique(sub.dmat@cols))
    cfeas <- length(unique(sub.dmat@rows))
  }

  max_dist <- suppressWarnings(max(sub.dmat@.Data))

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
  epsilon_lower_lim <- maxdist / (.Machine$integer.max/64 -2)
  epsilon <- if (tolerance>0 & rfeas>1 & cfeas>1) {
    max(epsilon_lower_lim, tolerance/(rfeas + cfeas - 2))
  } else epsilon_lower_lim
  options(old.o)
  return(epsilon)
}

parse_subproblems <- function(problems, min.controls,
                              max.controls, omit.fraction, hints,
                              total.n, TOL.in, resolution = NULL)
{
  np <- length(problems)
  subproblemids <- names(problems)
  if (is.null(subproblemids)) subproblemids  <- character(1L)
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

  #when no hint is provided, hints is a list of length np where everything is null
  #when hint is provided, hint is a list of NodeInfo elements
  hint.provided <- all(sapply(hints, function(x) !is.null(x)))

  hint.list <- mapply(prepare_subproblem_hint,
                      d = problems,
                      hint = hints,
                      SIMPLIFY = FALSE)

  if (!is.null(resolution))
  {
    length.res <- length(resolution)
    names.res <- names(resolution)

    if (all(is.null(names.res)))
    {
      if (identical(np, 1L) &&
          identical(length.res, 1L))
      {
        names.res <- character(1L)
        names(resolution) <- names.res
      } else {
        stop("No subproblem ids are specified.")
      }
    }

    if (hint.provided)
    {
      if(any(!names.res %in% subproblemids))
      {
        stop("Group specified in resolution argument cannot be found.")
      }
      # if resolution is specified for a few subproblems, do the default behavior first:
      tol.fracs <- mapply(get_tol_frac,
                          d = problems,
                          total_n = rep(total.n, np),
                          np = rep(np, np))
      tolerances <- tol.fracs * TOL.in

      epsilons <- mapply(get_epsilon,
                         hint = hint.list,
                         tolerance = tolerances,
                         subproblem = problems,
                         flipped = flippeds,
                         SIMPLIFY = FALSE)
      #Then, update with the new specified values
      epsilons[names.res] <- resolution

    } else {
      if(length.res != np ||
         !identical(sort(names.res), sort(subproblemids)))
      {
        stop("Resolutions must provided for each subproblem as a named list. Please ensure all subproblemids are correct and specified.")
      }


      if (length(resolution) == 1 &&
          names(resolution) == "" &&
          np == 1 && all(subproblemids == "")) # case of one subproblem
      {
        epsilons <- resolution
        names(epsilons) <- character(1L)
      } else {
        epsilons <- lapply(subproblemids, function(x) resolution[[x]])
        names(epsilons) <- subproblemids
      }

    }


  } else { #no epsilons provided, just use the default behavior and convert from tolerance
    tol.fracs <- mapply(get_tol_frac,
                        d = problems,
                        total_n = rep(total.n, np),
                        np = rep(np, np))
    tolerances <- tol.fracs * TOL.in

    epsilons <- mapply(get_epsilon,
                       hint = hint.list,
                       tolerance = tolerances,
                       subproblem = problems,
                       flipped = flippeds,
                       SIMPLIFY = FALSE)
  }

  tmp <- list(subproblems = problems,
              flipped_status = flippeds,
              min.controls = mincs,
              max.controls = maxcs,
              epsilons = epsilons,
              hints = hint.list)
  return(tmp)
}


find_subproblem_epsilons <- function(tol = NULL,
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


#' Get information about subproblems
#'
#' @param object optmatch object
#' @param type either "subproblem" or any of the column names from the subproblems table associated with an MCFSolutions object (e.g. "resolution")
#'
#' @return if \code{type = "subproblem"}, a named vector mapping units to subproblems is returned. If one or more columns from the subproblems table is specified, that column subset is returned. This might result in either a vector or data frame being returned.
#' @export
#'
#' @examples
#' data(nuclearplants)
#' f1 <- fullmatch(pr ~ cost + strata(pt), data=nuclearplants)
#' get_subproblem_info(object = f1, type = c("resolution", "groups"))
#' get_subproblem_info(object = f1, type = "subproblem")
get_subproblem_info <- function(object, type)
{
  if (!inherits(object, "optmatch")){
    stop("Please provide an optmatch object")
  }
  if (all(type == "subproblem"))
  {
    return(attr(object, "subproblem"))
  } else {
    subprob.attributes <- colnames(attr(object, "MCFSolutions")@subproblems)

    idx <- pmatch(type, subprob.attributes, nomatch = NA)
    idx <- idx[!is.na(idx)]

    tb <- attr(object, "MCFSolutions")@subproblems
    res <- tb[, idx]
    return(res)
  }
}

#' Parse hints for metadata
#'
#' @param full.hint
#'
#' @return
#' @export
#'
#' @examples
parse_hint_metadata <- function(full.hint)
{
  # assuming that the hint is not null
  # also assuming existence of MCFSolutions object

  groupids <- get_subproblem_info(full.hint, type = "groups")
  feasibility <- get_subproblem_info(full.hint, type = "feasible")
  names(feasibility) <- groupids
  resolutions <- get_subproblem_info(full.hint, type = "resolution")
  names(resolutions) <- groupids
  regret_gap <- get_subproblem_info(full.hint, type = "primal_value") - get_subproblem_info(full.hint, type = "dual_value")
  names(regret_gap) <- groupids

  rl <- list(feasible = feasibility,
       resolution = resolutions,
       regret = regret_gap)
  return(rl)
}
