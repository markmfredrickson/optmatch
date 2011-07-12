maxControlsCap <- function(distance, min.controls=NULL, subclass.indices=NULL)
{
  distance <- as.matrix(distance) # turns an ISM into a matrix, temporary cast
############################################################
# CHECK DIMNAMES OF DISTANCE			   #
############################################################
if (!is.list(distance) & !is.matrix(distance))
  stop("argument \'distance\' must be a matrix or list") 
if (is.matrix(distance))
  {
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
      idc <- factor(rep("m", length(nmtrt)+length(nmctl)))
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
mncpt <- fullmatchNumControlsHandling(min.controls, levels(idc), "min.controls")
############################################################
# HANDLE OMIT.FRACTION	                                   #
############################################################
omf <- rep(NA, nlevels(idc))
names(omf) <- levels(idc)

############################################################
# CREATE STRATUM IDENTIFIER				   #
############################################################
#ssizes <- lapply(split(idc, idc), length)

################################################################
# MARK AND SEPARATE UNITS BELONGING TO STRATA LACKING ROW UNITS#
# OR TO STRATA LACKING CONTROL UNITS			       #
################################################################
inrow <- (rns %in% nmtrt)
incol <- (rns %in% nmctl)

nctls <- unsplit(tapply(incol,idc,
                        function(x){rep(sum(x),length(x))}), 
                 idc)
ntrs <-  unsplit(tapply(inrow,idc,
                        function(x){rep(sum(x),length(x))}), 
                 idc)
mgrp <- (nctls>0) & (ntrs>0)

rnl <- split(rns[(mgrp&inrow)], idc[(mgrp&inrow)])
cnl <- split(rns[(mgrp&incol)], idc[(mgrp&incol)])
sfs <- unique(as.character(idc[mgrp]))

lmxc <- rep(Inf, nlevels(idc))
names(lmxc) <- levels(idc)

if (is.null(min.controls))
   {
   gmnc <- rep(0, nlevels(idc))
   names(gmnc) <- levels(idc)
   } else
   {
   gmnc <- mncpt
   gmnc[!(names(gmnc) %in% sfs)] <- 0
   }

for (i in sfs)
    {
    trnl <- rnl[[i]]
    tcnl <- cnl[[i]]
    tlmxc <- lmxc[i]
    tgmnc <- gmnc[i]
    if (is.matrix(distance))
      {
      tdm <- (distance[trnl,tcnl]>0)/(distance[trnl,tcnl]<Inf)
      } else
    {tdm <- (distance[[i]][trnl,tcnl]>0)/(distance[[i]][trnl,tcnl]<Inf)}
    tdm <- matrix(tdm, length(trnl), length(tcnl),
                  dimnames=list(trnl, tcnl))
# FEASIBILITY CHECK -- temp depends on whether problem requires flipping
    if (switch(1+is.na(omf[i]), omf[i]>0,  TRUE ))
      {
        ncol <- length(tcnl)
        nrow <- length(trnl)
        temp <- SubDivStrat(rownames=trnl, colnames=tcnl, distmat=tdm,
       max.cpt=min(tlmxc, ncol),
       min.cpt=max(tgmnc, 1/nrow), tolerance=.5, 
       omit.fraction=switch(1+is.na(omf[i]), omf[i], NULL))
      } else
    {
        ncol <- length(trnl)
        nrow <- length(tcnl)
        temp <- SubDivStrat(rownames=tcnl, colnames=trnl, distmat=t(tdm),
       max.cpt=min(1/tgmnc, ncol),
       min.cpt=max(1/tlmxc, 1/nrow), tolerance=.5, 
       omit.fraction=switch(1+is.na(omf[i]), -omf[i], NULL))

    }
# IF THE PROBLEM IS FEASIBLE, SET TLMXC TO GREATEST OBTAINED
# RATIO OF CONTROLS TO TREATED UNITS.  THIS MAY BE MUCH LESS THAN
# THE GENERIC BOUND WE'D OTHERWISE USE.
    if (!all(is.na(temp$cells))&& !all(temp$cells=="NA") && !all(temp$cells=="0"))
       {
       tlmxc <- max(apply(
                      table(temp$cells[temp$cells!='0'],
                            names(temp$cells)[temp$cells!='0'] %in% trnl),
                        1, function(x) {x[1]/x[2]}), na.rm=TRUE)
     }
    
    # ONLY GO FURTHER IF LEAST RESTRICTIVE TLMXC GAVE FEASIBILITY
    # ALSO, FOR THE TIME BEING, NEGATIVE OMIT.FRACTION NOT DEALT WITH
    if (!all(is.na(temp$cells))&& !all(temp$cells=="NA") && !all(temp$cells=="0") &&
        switch(1+is.na(omf[i]), omf[i]>=0, TRUE))
      {
    
    if (tgmnc < 1)
      {
        # SHOULD TLMXC ALSO BE SET TO ONE OR LESS?
        ncol <- length(trnl)
        nrow <- length(tcnl)
        temp <- SubDivStrat(rownames=tcnl, colnames=trnl, distmat=t(tdm),
                            max.cpt=min(1/tgmnc, ncol), min.cpt=1,
                            tolerance=.5, omit.fraction=
                            switch(1+is.na(omf[i]), -omf[i], NULL))
        flipflag <- !all(is.na(temp$cells)) && !all(temp$cells=="NA") && !all(temp$cells=="0")
      } else {flipflag <- FALSE}
 
    if (flipflag)
      {
      # CASE THAT BOTH GMNC AND OPTIMUM LMXC ARE LESS THAN 1
      # (BUT FIRST, CHECK THAT TLMXC ISN'T ALREADY AS LARGE AS IT CAN GET)
        if (min(1/tgmnc,length(trnl))!=max(1, 1/tlmxc))
        {
        tlmxc <- 
       optimize( function(invlmxc, rown1, coln1, dist1, gmnc1, omf1) {
       ifelse(!all(SubDivStrat(rownames=coln1, colnames=rown1, distmat=t(dist1),
       max.cpt=min(1/gmnc1, length(rown1)), min.cpt=invlmxc,
       tolerance=.5, omit.fraction= switch(1+is.na(omf[i]), -omf[i],
       NULL) )$cells=="NA") ,
              invlmxc, -invlmxc)
                  },
                upper=min(1/tgmnc,length(trnl)), lower=max(1, 1/tlmxc),
                tol=1, maximum=TRUE, rown1=trnl, coln1=tcnl,
                dist1=tdm, gmnc1=tgmnc, omf1=omf[i]
                )$objective
        tlmxc <- 1/floor(tlmxc)
      } 
      } else
       {  
       # TREAT USUAL CASE, LMXC WILL BE SOMEWHERE ABOVE MAX(GMNC, 1)
      # (BUT FIRST, CHECK THAT TLMXC ISN'T ALREADY AS LARGE AS IT CAN GET)
       if (max(tgmnc,1)!=min(length(tcnl), tlmxc))
         {
       tlmxc <- ceiling(
       optimize( function(lmxc1, rown1, coln1, dist1, gmnc1, omf1) {
       ifelse(!all(SubDivStrat( rownames=rown1, colnames=coln1, distmat=dist1,
       min.cpt=max(gmnc1, 1/length(rown1)), max.cpt=lmxc1,
       tolerance=.5, omit.fraction= switch(1+is.na(omf[i]), omf[i],
       NULL) )$cells=="NA") ,
              lmxc1, 2*length(coln1) - lmxc1)
                  },
                lower=max(tgmnc,1), upper=min(length(tcnl), tlmxc), tol=1,
                maximum=FALSE, rown1=trnl, coln1=tcnl,
                dist1=tdm, gmnc1=tgmnc, omf1=omf[i]
                )$objective )
       } }
  lmxc[i] <- tlmxc
  }
  }
list(given.min.controls=gmnc, strictest.feasible.max.controls=lmxc)

}
 
