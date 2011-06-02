makedist <- function(structure.fmla, data,
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
