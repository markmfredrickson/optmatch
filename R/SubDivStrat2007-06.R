SubDivStrat <- function(rownames, colnames, distmat, min.cpt,
max.cpt, tolerance, omit.fraction=NULL, matched.distances=FALSE)
{
if (min.cpt <=0 | max.cpt<=0) stop("inputs min.cpt, max.cpt must be positive")

if (!all(rownames %in% dimnames(distmat)[[1]]))
   stop("input \'rownames\' may only contain row names of input \'distmat\'")
if (!all(colnames %in% dimnames(distmat)[[2]]))
   stop("input \'rownames\' may only contain col. names of input \'distmat\'")

rownames <- as.character(rownames)
colnames <- as.character(colnames)
if (length(rownames)>1)
  {
  dm <- as.matrix(distmat[rownames, colnames])
  } else
  {
  dm <- matrix(distmat[rownames, colnames], length(rownames),
               length(colnames), dimnames=list(rownames, colnames))
  }
rfeas <- apply(dm,1,function(x) {any(is.finite(x))})
cfeas <- apply(dm,2,function(x) {any(is.finite(x))})
if (floor(min.cpt) > ceiling(max.cpt) | ceiling(1/min.cpt) < floor(1/max.cpt)) 
   {
   ans <- rep("NA",length(rownames)+length(colnames))
   names(ans) <- c(rownames, colnames)
   return(list(cells=ans, maxerr=NULL, distance=NULL))
   }

if (is.null(omit.fraction))
   {
   f.ctls <- 1
   dm <- matrix(dm[rfeas, cfeas], sum(rfeas), sum(cfeas),
	 dimnames=list(rownames[rfeas], colnames[cfeas]))
   } else 
   { 
   if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction
      > 1) {stop("omit.fraction must be null or between 0 and 1")}
   f.ctls <- 1-omit.fraction 
   dm <- matrix(dm[rfeas,], sum(rfeas), length(colnames),
	 dimnames=list(rownames[rfeas], colnames))
   }

if (any(rfeas) & any(cfeas))
   {
   old.o <- options(warn=-1)
   if (any(dm > 0 & is.finite(dm))) {
   reso <- (.Machine$integer.max/64 -2)/max(dm[is.finite(dm)]) } else {
	reso <- min(.Machine$integer.max/64 -2, 
		    (sum(rfeas)+sum(cfeas))/tolerance)   }
   if (tolerance>0 & sum(rfeas)>1 & sum(cfeas)>1) 
      {
   reso <- min(reso, (sum(rfeas) + sum(cfeas) - 2)/tolerance) 
      }
   options(old.o)
   temp <- fmatch(floor(dm*reso), max.row.units=ceiling(1/min.cpt), 
		  max.col.units=ceiling(max.cpt),
		min.col.units=max(1, floor(min.cpt)), f=f.ctls) 
   if (any(temp<0)) {maxerr <- 0} else 
      {
      maxerr <- sum(temp*dm, na.rm=TRUE) - 
	     sum(temp*floor(dm*reso), na.rm=TRUE)/reso +
	     (sum(rfeas)>1 & sum(cfeas)>1)*
	     (sum(rfeas) + sum(cfeas) - 2 - sum(temp))/reso 
      }
if ( maxerr > tolerance)
   {
   temp1 <- temp
   temp2 <- fmatch(round(dm*reso), max.row.units=ceiling(1/min.cpt), 
		  max.col.units=ceiling(max.cpt),
		min.col.units=max(1, floor(min.cpt)), f=f.ctls) 
   if  (sum(temp1*dm, na.rm=TRUE)<=sum(temp2*dm, na.rm=TRUE)) 
        {temp <- temp1} else {temp <- temp2}
     maxerr <- sum(temp*dm, na.rm=TRUE) - sum(temp1*floor(dm*reso), na.rm=TRUE)/reso +
	     (max(1, sum(rfeas)-1)+ max(1, sum(cfeas)-1) - 
	      (sum(rfeas)==1 & sum(cfeas)==1) - sum(temp1))/reso 
    }
if (matched.distances)
  {
if (all(temp>=0))
   {
   dma <- max(dm[as.logical(temp)])
   dist <- c(apply(temp*pmin(dm, dma), 1, sum), 
		apply(temp*pmin(dm, dma), 2, sum))
dist[c(rep(FALSE, dim(temp)[1]), 
       apply(temp*apply(temp, 1, sum), 2, sum)==1)] <- NA
   dist[c(apply(temp, 1, sum)>1, apply(temp, 2, sum)>1)] <- NA
   } else
   {
   dist <- rep(NA, sum(dim(temp)))
   mode(dist) <- "numeric"
   names(dist) <- unlist(dimnames(temp))
   }
 } else {dist <- 0}
   } else { temp <- 0 ; maxerr <- 0 ; dist <- 0}
stratumlabels <- function(matches)
{
rdm <- dim(matches)[1]
cdm <- dim(matches)[2]
## LABEL CONTROL SUBJECTS ACCORDING TO ROW NUMBER OF FIRST TREATED 
## SUBJECT TO WHOM THEY ARE MATCHED
clbl <- apply(matches, 2, function(x)
                      { (c(1:rdm, 0)[c(as.logical(x), TRUE)])[1] } )

c(apply(rep(clbl, rep(rdm, cdm)) * matches, 1, function(x) {
  c(x, 0)[c(as.logical(x), TRUE)][1]} ), clbl)

}

if (all(temp<=0) & any(temp<0))
{
ans <- rep("NA",length(rownames)+length(colnames))
names(ans) <- c(rownames, colnames)
} else
{
ans <- rep("0",length(rownames)+length(colnames))
names(ans) <- c(rownames, colnames)
if (any(temp>0))
   {
   dimnames(temp) <- dimnames(dm)
   mout <- stratumlabels(temp)
   ans[names(mout)] <- mout
   }
}
#if (any(ans=="0")) {ans[ans=="0"] <- paste("0", 1:sum(ans=="0"), sep="")}
if (any(ans=="0")) {ans[ans=="0"] <- NA}

list(cells=ans, err=maxerr, match.distance=dist)
}
