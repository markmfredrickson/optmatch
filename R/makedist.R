################################################################################
# makedist: produces DistanceSpecificaitons from data and a distance function
################################################################################

# an internal function to handle the heavy lifting
# z is a treatment indicator for the data
#
makedist <- function(z, data, distancefn, within = NULL) {
  # conversion and error checking handled by toZ
  z <- toZ(z)

  # next, data should be same length/size as z
  # I would have liked to have "duck typed" data by using
  # hasMethod(f = "dim", signature = class(data))
  # but "vector" has a dim method (though it returns null!), which causes this to fail
  if (inherits(data, "matrix") | inherits(data, "data.frame")) {
    data.len <- dim(data)[1]
  } else {
    data.len <- length(data)
  }

  if (data.len != length(z)) {
    stop("Treatment indicator and data must have the same length.")
  }

  # data can be either a vector or a data.frame/matrix
  if (is.vector(data)) {
    namefn <- names
  } else {
    namefn <- rownames
  }
  rns <- namefn(data)[as.logical(z)]
  cns <- namefn(data)[!as.logical(z)]

  # rns <- namefn(subset(data, as.logical(z))) # the d.f subset requries the as.logical. weird
  # cns <- namefn(subset(data, !as.logical(z)))

  if (length(cns) == 0 | length(rns) == 0) {
    stop(paste("Data must have ", ifelse(is.vector(data), "names", "rownames"), ".", sep = ""))
  }

  warning.requested <- getOption("optmatch_warn_on_big_problem", default = TRUE)

  if (is.null(within)) {
    # without a within, make a dense matrix
    nc <- length(cns)
    nr <- length(rns)
    # matrices have column major order
    treatmentids <- rep(rns, nc)
    controlids <- rep(cns, each = nr)

    if ((as.numeric(nc) * nr > getMaxProblemSize()) && warning.requested) {

      warning("match_on has been asked to compute a large number of treatment-control
distances.  You can reduce this number by providing an appropriate 'within'
argument; see 'match_on', 'exactMatch' and 'caliper' documentation for
details.  Alternatively you could increase the problem size limit; see
documentation of 'getMaxProblemSize' and 'setMaxProblemSize'.")

    }

  } else {
    # with a within, make a copy and only fill in the finite entries of within

    if (!(all(rns %in% within@rownames)) | !(all(cns %in% within@colnames))) {
      stop("Row and column names of within must match those of the data.")
    }

    if(!all(within@rownames %in% rns) | !all(within@colnames %in% cns))
    {
        if(!is(within, "BlockedInfinitySparseMatrix")) {
          within <- subset(within, within@rownames %in% rns, within@colnames %in% cns)
        } else {
          tmpISM <- subset(within, within@rownames %in% rns, within@colnames %in% cns)
          tmpBISM <- as(tmpISM, "BlockedInfinitySparseMatrix")
          newgroups <- within@groups[names(within@groups) %in% c(tmpBISM@rownames, tmpBISM@colnames)]
          tmpBISM@groups <- newgroups
          names(tmpBISM@groups) <- names(newgroups)
          within <- tmpBISM
        }

    }

    treatmentids <- within@rownames[within@rows]
    controlids <- within@colnames[within@cols]



    # TODO: check that the rownames, colnames of within match data
    subprobs <- findSubproblems(within)
    issue.warning <- FALSE # just issue 1 warning, even if multiple subs fail
    for (s in subprobs) {
      issue.warning <- issue.warning || (dim(prepareMatching(s))[1] > getMaxProblemSize())
    }

    if(issue.warning && warning.requested) {
      warning("I've been asked to compute a large number of treatment-control distances.
Even with the provided 'within' argument, the result will present too large
an optimization problem for optimal matching. Please use a more restrictive
'within' specification, either in a revision of this call or a follow-up
elaborating on it; see 'exactMatch' and 'caliper' documentation for details.")
    }
  }

    dists <- distancefn(cbind(treatmentids, controlids), data, z)
    if (length(dists)!=length(treatmentids))
        stop("Internal optmatch error - malformed distances.\n (Try a different method argument to match_on.)")

  # z was copied <- toZ(z) so should be safe to rm
  # before massive matrix alloc
  rm(treatmentids, controlids, z)

  if(is.null(within)) {
      res <- new("DenseMatrix", matrix(dists, nrow = nr, ncol = nc, dimnames =
                                       list(treatment = rns, control = cns)))
  } else {
    # using `replace` will cause issues down the road with `[<-`
    res <- within
    res@.Data <- dists
  }
  return(res)
}


makedistold <- function(structure.fmla, data,
                     fn=function(trtvar, dat, ...){
                       matrix(0, sum(trtvar), sum(!trtvar),
                              dimnames=list(names(trtvar)[trtvar],
                                names(trtvar)[!trtvar]))},
                     ...)
  {
  if (!attr(terms(structure.fmla), "response")>0)
    stop("structure.fmla must specify a treatment group variable")
  fn <- match.fun(fn)

### WHEN THIS FUNCTION IS WRAPPED TO, THIS IS HOW INFO ABOUT WHICH
### GENERATION PARENT FRAME structure.fmla IS TO BE EVALUATED IN IS PASSED
  pframe.generation <- 1
  if (!is.null(attr(structure.fmla, "generation.increment")))
    pframe.generation <- pframe.generation +
      attr(structure.fmla, "generation.increment")

  zpos <- attr(terms(structure.fmla), "response")
  vars <- eval(attr(terms(structure.fmla), "variables"), data,
             parent.frame(n=pframe.generation))
  zzz <- vars[[zpos]]
  if (!is.numeric(zzz) & !is.logical(zzz))
    stop("treatment variable (LHS of structure.fmla) must be numeric or logical")
  if (any(is.na(zzz)))
      stop("NAs not allowed in treatment variable (LHS of structure.fmla)")
  if (all(zzz>0))
    stop("there are no controls (LHS of structure.fmla >0)")
  if (all(zzz<=0))
    stop("there are no treatment group members (LHS of structure.fmla <=0)")

  zzz <- (zzz>0)
  vars <- vars[-zpos]
  names(zzz) <- row.names(data)
if (length(vars)>0)
  {
  ss <- interaction(vars, drop=TRUE)
} else ss <- factor(zzz>=0, labels="m")
  ans <- tapply(zzz, ss, FUN=fn,
                dat=data, ..., simplify=FALSE)
  FUNchk <- unlist(lapply(ans,
                          function(x){!is.matrix(x) & !is.vector(x)}))
  if (any(FUNchk)) { stop("fn should always return matrices")}

  mdms <- split(zzz,ss)
  NMFLG <- FALSE

  for (ii in (seq_along(ans)))
  {
    dn1 <- names(mdms[[ii]])[mdms[[ii]]]
    dn2 <- names(mdms[[ii]])[!mdms[[ii]]]
    if (is.null(dim(ans[[ii]])))
      {
      if (length(dn1)>1 & length(dn2)>1)
        { stop("fn should always return matrices")}
      if (length(ans[[ii]])!=max(length(dn1), length(dn2)))
        { stop(paste("unuseable fn value for stratum", names(ans)[ii]))}

      if (is.null(names(ans[[ii]])))
          {
           ans[[ii]] <-  matrix(ans[[ii]], length(dn1), length(dn2),
                         dimnames=list(dn1,dn2))
         } else
        { if (length(dn1)>1)
            {
              ans[[ii]] <- matrix(ans[[ii]],length(dn1), 1,
                                 dimnames=list(names(ans[[ii]]), dn2))
            } else {
              ans[[ii]] <- matrix(ans[[ii]], 1, length(dn2),
                                 dimnames=list(dn1, names(ans[[ii]])))
            }
        }
      } else {
      if (!all(dim(ans[[ii]])==c(length(dn1), length(dn2))))
        { stop(paste("fn value has incorrect dimension at stratum",
                     names(ans)[ii])) }
      if (is.null(dimnames(ans[[ii]])))
        {
          dimnames(ans[[ii]]) <- list(dn1, dn2)
          NMFLG <- TRUE
        } else {
          if (!all(dn1%in%dimnames(ans[[ii]])[[1]]) |
              !all(dn2%in%dimnames(ans[[ii]])[[2]]) )
            { stop(paste(
                    "dimnames of fn value don't match unit names in stratum",
                         names(ans)[ii])) }
        }

      }
  }

  if (NMFLG){
warning("fn value not given dimnames; assuming they are list(names(trtvar)[trtvar], names(trtvar)[!trtvar])")}


  attr(ans, 'row.names') <- names(zzz)
  class(ans) <- c('optmatch.dlist', 'list')
  ans
  }
