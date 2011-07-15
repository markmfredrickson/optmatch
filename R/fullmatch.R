fullmatch <- function(distance,  
    min.controls = 0, 
    max.controls = Inf, 
    omit.fraction = NULL, 
    tol = .001) {

  ### Checking Input ###
  
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification object")      
  }

  # we expect the following to be defined for the distance object
  # put any functions in this list that are called directly on distance
  methods <- c("dim", "dimnames", "prepareMatching", "subproblems",
    "is.numeric")
  dist.class <- class(distance)
  lapply(methods, function(m) {
    if (!hasMethod(m, dist.class)) {
      # skip the FUN = ... in the call stack
      stop(paste("argument \'distance\' must have a", m, "method."), call. = F)  
    }
  })

  dnms <- dimnames(distance)
  if (is.null(dnms) | is.null(dnms[[1]]) | is.null(dnms[[2]])) {
    stop("argument \'distance\' must have dimnames") 
  }
  
  if (any(duplicated(unlist(dnms)))){ 
    stop("dimnames of argument \'distance\' contain duplicates") 
  }

  nmtrt <- dnms[[1]]
  nmctl <- dnms[[2]]

  # note: this next _should_ be unnecessary, the objects should do this
  # but better safe than sorry
  if (!identical(dim(distance), c(length(nmtrt), length(nmctl)))) {
    stop("argument \'distance\' dimensions do not match row and column names")  
  }

  if (!is.numeric(distance)) {
    stop("argument \'distance\' must be numeric")  
  }

  # problems is guaranteed to be a list of DistanceSpecifictions
  # it may only have 1 entry
  problems <- findSubproblems(distance)

  # the number of problems should match the argument lengths for
  # min, max, and omit

  np <- length(problems)
  if (length(min.controls) > 1 & np != length(min.controls)) {
    stop(paste("Length of \'min.controls\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  if (length(max.controls) > 1 & np != length(max.controls)) {
    stop(paste("Length of \'max.controls\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  if (!is.null(omit.fraction) & length(omit.fraction) > 1 & np !=
    length(omit.fraction)) {
    stop(paste("Length of \'omit.fraction\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  
  # reset the arguments to be the right length if they are not
  if (length(min.controls) == 1) {
    min.controls <- rep(min.controls, np)   
  }
  if (length(max.controls) == 1) {
    max.controls <- rep(max.controls, np)   
  }
  if (length(omit.fraction) == 1) {
    omit.fraction <- rep(omit.fraction, np)   
  }

  idc <- factor(rep("m", length(nmtrt) + length(nmctl)))
  if (any(is.na(suppressWarnings(as.numeric(c(nmtrt,nmctl)))))) {
    names(idc) <- sort(c(nmtrt,nmctl))
  } else {
    names(idc) <- c(nmtrt,nmctl)[order(as.numeric(c(nmtrt,nmctl)))]
  }
  rns <- names(idc)

############################################################
# HANDLE DIFFERENT INPUT FORMS FOR MIN.CONTROLS		   #
############################################################
  mncpt <- fullmatchNumControlsHandling(min.controls,levels(idc),"min.controls")
############################################################
# HANDLE DIFFERENT INPUT FORMS FOR MAX.CONTROLS		   #
############################################################
  mxcpt <- fullmatchNumControlsHandling(max.controls,levels(idc),"max.controls")
############################################################
# HANDLE DIFFERENT INPUT FORMS FOR OMIT.FRACTION	   #
############################################################
  if (is.null(omit.fraction)) 
  {
    omf <- rep(NA, nlevels(idc))
      names(omf) <- levels(idc)
  } else
  {
    if (any(abs(omit.fraction)>1, na.rm=TRUE) | !is.numeric(omit.fraction))
    { stop("omit.fraction must be NULL or numeric between -1 and 1") }
    if (length(omit.fraction)==1)
    {
      omf <- rep(omit.fraction, nlevels(idc)) 
        names(omf) <- levels(idc)
    } else
    {
      if (!all(levels(idc) %in% names(omit.fraction))) {
        stop("\'omit.fraction\' not specified for some subclasses") }
      omf <- omit.fraction[levels(idc)]
    } 
    if (any(omf > 0 & mxcpt <= .5, na.rm=TRUE) )
      stop("positive \'omit.fraction\' with \'max.controls\' <= 1/2 not permitted")
        if (any(omf < 0 & mncpt >= 2, na.rm=TRUE) )
          stop("negative \'omit.fraction\' with \'min.controls\' >= 2 not permitted")
  }
############################################################
# CREATE STRATUM IDENTIFIER				   #
############################################################
#strat <- character(length(rns))
##strat <- as.character(idc) # NEW
#names(strat) <- rns
#ssizes <- lapply(split(strat, idc), length)
#split(strat, idc) <- rep(names(ssizes), unlist(ssizes))


################################################################
# MARK AND SEPARATE UNITS BELONGING TO STRATA LACKING ROW UNITS#
# OR TO STRATA LACKING CONTROL UNITS			       #
################################################################
  inrow <- (rns %in% nmtrt)
  incol <- (rns %in% nmctl)

#nctls <- ntrs <- numeric(length(idc))
#split(nctls, idc) <- rep(unlist(lapply(split( 
#			   incol, idc),sum)), unlist(ssizes))
#split(ntrs, idc) <- rep(unlist(lapply(split(
#			   inrow, idc), sum)), unlist(ssizes))
  nctls <- unsplit(tapply(incol,idc,
        function(x){rep(sum(x),length(x))}), 
      idc)
  ntrs <-  unsplit(tapply(inrow,idc,
        function(x){rep(sum(x),length(x))}), 
      idc)
  mgrp <- (nctls>0) & (ntrs>0)
  strat.abv <- abbreviate(as.character(idc), 2)
  names(strat.abv) <- names(idc)
#strat.abv[!mgrp] <- paste(strat.abv[!mgrp], "0", sep=".")
  strat.abv[!mgrp] <- NA

  rnl <- split(rns[(mgrp&inrow)], factor(idc[(mgrp&inrow)]))
  cnl <- split(rns[(mgrp&incol)], factor(idc[(mgrp&incol)]))
  sfs <- levels(factor(idc[mgrp]))
  err <- numeric(length(sfs)) ; names(err) <- sfs
  TOL <- tol*sum(mgrp)

# TODO: remove cast when removing all list capabilities
  distance <- as.matrix(distance)

  for (i in sfs)
  { 
    if (switch(1+is.na(omf[i]), omf[i]>0,  mxcpt[i] > .5 ))
    {
      nrow <- length(rnl[[i]])
        ncol <- length(cnl[[i]])
        tol.frac <- (nrow+ncol-2)/(sum(mgrp)-2*length(rnl))
        temp <- SubDivStrat(rownames=rnl[[i]], colnames=cnl[[i]], 
            distmat=switch(mode(distance),list=distance[[i]],numeric=distance),
            max.cpt=min(mxcpt[i], ncol), 
            min.cpt=max(mncpt[i], 1/nrow), tolerance=(TOL*tol.frac), 
            omit.fraction=switch(1+is.na(omf[i]), omf[i], NULL))
    } else
    {
      ncol <- length(rnl[[i]])
        nrow <- length(cnl[[i]])
        tol.frac <- (nrow+ncol-2)/(sum(mgrp)-2*length(rnl))
        temp <- SubDivStrat(rownames=cnl[[i]], colnames=rnl[[i]], 
            distmat=t(switch(mode(distance),list=distance[[i]],numeric=distance)),
            max.cpt=min(1/mncpt[i], ncol), 
            min.cpt=max(1/mxcpt[i], 1/nrow), tolerance=(TOL*tol.frac), 
            omit.fraction=switch(1+is.na(omf[i]), -omf[i], NULL))
    }

    strat.abv[names(temp$cells)] <-
      ifelse(is.na(temp$cells),NA,
          paste(strat.abv[names(temp$cells)],
            temp$cells, sep=".") )

      if (!any(!is.na(temp$cells) & temp$cells=="NA")) 
      {
        err[i] <- temp$err
      }
    NULL
  }
    strat.abv <- as.factor(strat.abv)
      if (inherits(distance, "optmatch.dlist"))
      {
        if (all(attr(distance, "row.names")%in%names(strat.abv)))
        {
          inrow <- inrow[match(attr(distance, "row.names"), names(strat.abv))]
            strat.abv <- strat.abv[match(attr(distance, "row.names"), names(strat.abv))]
        } else {
          warning("row.names attribute of distance doesn't match dimnames of dist matrices")
        }
      }

    class(strat.abv) <- c("optmatch", "factor")
      attr(strat.abv, "exceedances") <- err
      if (sum(err, na.rm=TRUE)>TOL) 
      {
        warning(
            paste("prescribed tol of ", tol, "per obs. poss. exceeded by up to ", 
              round(sum(err), 3), ".", sep="") )
      }
    attr(strat.abv, "call") <- match.call()

      attr(strat.abv, "contrast.group") <- inrow
      attr(strat.abv, "matched.distances") <- matched.distances(strat.abv, distance)
      strat.abv
}
