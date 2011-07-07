fullmatch <- function(distance,  min.controls=0, max.controls=Inf, 
omit.fraction=NULL, tol=.001, subclass.indices=NULL)
{
############################################################
# CHECK DIMNAMES OF DISTANCE			   #
############################################################
if (!is(distance, "DistanceSpecification"))
  stop("distance must be a DistanceSpecification object")

# TODO: remove cast when removing all list capabilities
distance <- as.matrix(distance)

if (is.matrix(distance))
  {
    if (is.logical(distance))
      distance <- matrix(as.numeric(distance),
                         nrow(distance), ncol(distance),
                         dimnames=dimnames(distance))
    
  if (!is.numeric(distance))
    stop("matrix \'distance\' must be of mode numeric")
  if (is.null(dimnames(distance))) {
   stop("argument \'distance\' must have dimnames") }
  if (any(duplicated(unlist(dimnames(distance)))))
   { stop("dimnames of argument \'distance\' contain duplicates") }
  nmtrt <- dimnames(distance)[[1]]
  nmctl <- dimnames(distance)[[2]]
} else
{
if (!all(unlist(lapply(distance,function(x){is.matrix(x) | is.null(x)}))))
   {
   bads <- (1:length(distance))[!sapply(distance,function(x){
                                                  is.matrix(x) | is.null(x)})]
  stop(paste("elements", 
             ifelse(length(bads)>1, paste(bads[1],"...",sep=""),bads[1]), 
             "of list \'distance\' fail to be (numeric) matrices.",
             sep=" "))
   }
distance[sapply(distance,is.null)] <- NULL
if (any(lcl <- sapply(distance,is.logical)))
  distance[lcl] <- lapply(distance[lcl],
                          function(x) matrix(as.numeric(x),
                                                nrow(x), ncol(x),
                                                dimnames=dimnames(x))
                          )

if (!all(unlist(lapply(distance,is.numeric))))
  stop("elements of list \'distance\' must be of mode numeric")
if (any(unlist(lapply(distance, function(x){is.null(dimnames(x))}))))
  stop("matrices in list \'distance\' must have dimnames")
  nmtrt <- unlist(lapply(distance, function(x){dimnames(x)[[1]]}))
  nmctl <- unlist(lapply(distance, function(x){dimnames(x)[[2]]}))
}
############################################################
# HANDLE DIFFERENT INPUT FORMS FOR SUBCLASS.INDICES	   #
############################################################
if (is.matrix(distance))
  {
if (is.data.frame(subclass.indices))
   {
   if (is.null(row.names(subclass.indices)) | 
      !all(row.names(subclass.indices) %in%
      c(nmtrt,nmctl)))
      stop("row names of data frame \'subclass.indices\' must exist and occur in dimnames of distance")
   rns <- row.names(subclass.indices)
   subclass.indices <- interaction(subclass.indices, drop=TRUE)
   names(subclass.indices) <- rns
   }
if (is.factor(subclass.indices)) 
   {
   if (is.null(names(subclass.indices)) | 
      !all(names(subclass.indices) %in% c(nmtrt,nmctl)) ) 
      stop("names of factor \'subclass.indices\' must exist and occur in dimnames of distance")
   idc <- subclass.indices[!is.na(subclass.indices)]
   rns <- names(idc) 
   } else 
   {
   if (is.null(subclass.indices)) 
      {
      idc <- factor(rep("m",length(nmtrt)+length(nmctl)))
      if
      (any(is.na(suppressWarnings(as.numeric(c(nmtrt,nmctl))))))
      {names(idc) <- sort(c(nmtrt,nmctl))}
      else
      {names(idc) <- c(nmtrt,nmctl)[
		     order(as.numeric(c(nmtrt,nmctl)))]
      }
      rns <- names(idc)
      } else    
      {
      stop("argument \'subclass.indices\' must be a factor") 
      }
   }
} else
{
if (!is.null(subclass.indices))
  warning("argument \'subclass.indices\' ignored when \'distance\' is a list")
if (is.null(names(distance)))
  {
dnm <- paste("m", 1:length(distance), "l", sep="")
} else dnm <- names(distance)
dnm[names(distance)==""] <-
  paste("m", 1:length(distance), "l", sep="")[names(distance)==""]
names(distance) <- dnm
idc <- factor(c(rep(dnm, unlist(lapply(distance,function(x){dim(x)[1]}))),
                rep(dnm, unlist(lapply(distance,function(x){dim(x)[2]})))))
names(idc) <- c(nmtrt,nmctl)
if (any(is.na(suppressWarnings(as.numeric(c(nmtrt,nmctl))))))
   {idc <- idc[order(c(nmtrt,nmctl))]}
   else
   {idc <- idc[order(as.numeric(c(nmtrt,nmctl)))]
      }
rns <- names(idc)
}
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
