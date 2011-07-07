################################################################################
# makedist: produces DistanceSpecificaitons from data and a distance function
################################################################################

# an internal function to handle the heavy lifting
# z is a treatment indicator for the data
# 
makedist <- function(z, data, distancefn, exclusions = NULL) {
  
  # first, check z for correctness
  if (any(is.na(z))) {
    stop("NAs not allowed in treatment indicator.")  
  }

  if (length(unique(z)) != 2) {
    stop(paste("Treatment indicator must have exactly 2 levels not", length(levels(z))))  
  }

  if (is.factor(z)) {
    z <- as.numeric(z) - 1 # make a 0/1 numeric vector  
  }
  z <- as.logical(z) # at this point, z is numeric 0/1, so is safe to convert to logical

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

  rns <- namefn(subset(data, as.logical(z))) # the d.f subset requries the as.logical. weird
  cns <- namefn(subset(data, !as.logical(z)))

  if (length(cns) == 0 | length(rns) == 0) {
    stop(paste("Data must have ", ifelse(is.vector(data), "names", "rownames"), ".", sep = ""))  
  }

  if (is.null(exclusions)) {
    # without a exclusions, make a dense matrix    
    nc <- length(cns)
    nr <- length(rns)
    
    res <- matrix(0, nrow = nr, ncol = nc, dimnames = list(treatment = rns, control = cns))
    
    # matrices have column major order
    treatmentids <- rep(rns, nc)
    controlids <- rep(cns, each = nr)

  } else {
    # with a exclusions, make a copy and only fill in the finite entries of exclusions
    res <- exclusions
    
    if (!all(exclusions@rownames %in% rns) | !(all(rns %in% exclusions@rownames)) |
        !all(exclusions@colnames %in% cns) | !(all(cns %in% exclusions@colnames))) {
      stop("Row and column names of exclusions must match those of the data.")  
    }

    treatmentids <- res@rownames[res@rows]
    controlids <- res@colnames[res@cols]

    # TODO: check that the rownames, colnames of exclusions match data
  }

  # it would be nice if we could abstract this, like with subset(data, z),
  # I am unaware of a function that will do that.

  if (is.vector(data)) {
    dists <- distancefn(data[treatmentids], data[controlids])    
  } else {
    dists <- distancefn(data[treatmentids,, drop = FALSE], 
                        data[controlids,, drop = FALSE])          
  }

  res <- replace(res, 1:length(res), dists)

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

  for (ii in (1:length(ans)))
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
