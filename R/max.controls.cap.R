##' @export
##' @rdname minmaxctlcap
maxControlsCap <- function(distance, min.controls = NULL)
{
  # check if it is valid distance specification,
  # if not through an error message explaining the issue
  validDistanceSpecification(distance, stopOnProblem = TRUE)

  # if we get this far, the distspec is valid, and findSubproblems
  # will generate a list of also valid distspecs.
  sps <- findSubproblems(distance)
  distance <- lapply(sps, as.matrix) # turns an ISM into a matrix, temporary cast

  nmtrt <- rownames(distance)
  nmctl <- colnames(distance)

  # notes: idc is a vector where each entry is which subclass the unit is in
  # the names of idc map back to the units
  #

  # deleted a lot of the stuff surrounding subclass.indices, saving this
  # for reference
  # idc <- factor(c(rep(dnm, unlist(lapply(distance,function(x){dim(x)[1]}))),
  #              rep(dnm, unlist(lapply(distance,function(x){dim(x)[2]})))))
  # names(idc) <- c(nmtrt,nmctl)
  # if (any(is.na(suppressWarnings(as.numeric(c(nmtrt,nmctl))))))
  #    {idc <- idc[order(c(nmtrt,nmctl))]}
  #    else
  #    {idc <- idc[order(as.numeric(c(nmtrt,nmctl)))]
  #       }
  # rns <- names(idc)
  rns <- c(nmtrt, nmctl)

  # min.controls is valid if it is a scalar or a vector with the same names as sps
  if(length(min.controls) > 1 &
     (!all(names(min.controls) %in% names(sps)) |
     !all(names(sps) %in% names(min.controls)))) {
    stop("Names of 'min.controls' must match the subproblems. See 'findSubproblems' and 'exactMatch'.")
  }

  # if min.controls is null we set it to 0 by default
  if(is.null(min.controls)) {
    min.controls <- 0
  }

  # make it into a well named vector of min.control values
  if(length(min.controls) == 1) {
    min.controls <- rep(min.controls, length(sps))
    names(min.controls) <- names(sps)
  }

  .maxControlsCap <- function(p, mc) {
    p <- as.matrix(p) # MMF: easier to upgrade the function by using a matrix
    omf <- NA
    trnl <- rownames(p)
    tcnl <- colnames(p)
    tlmxc <- Inf
    tgmnc <- mc

    # MMF: this check is probably not required anymore: the prepareMatching()
    # method should take care of things that are entirely unmatchable.
    tdm <- (p > 0) / (p < Inf)
    tdm <- matrix(tdm, length(trnl), length(tcnl), dimnames = list(trnl, tcnl))

    # FEASIBILITY CHECK -- temp depends on whether problem requires flipping
    ncol <- length(tcnl)
    nrow <- length(trnl)

    temp <- SubDivStrat(rownames = trnl, colnames = tcnl, distspec = tdm,
      max.cpt = min(tlmxc, ncol),
      min.cpt = max(tgmnc, 1/nrow), tolerance=.5,
      omit.fraction = NULL)

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
        switch(1+is.na(omf), omf>=0, TRUE))
    {

      if (tgmnc < 1)
      {
# SHOULD TLMXC ALSO BE SET TO ONE OR LESS?
        ncol <- length(trnl)
          nrow <- length(tcnl)
          temp <- SubDivStrat(rownames=tcnl, colnames=trnl, distspec=t(tdm),
              max.cpt=min(1/tgmnc, ncol), min.cpt=1,
              tolerance=.5, omit.fraction=
              switch(1+is.na(omf), -omf, NULL))
          flipflag <- !all(is.na(temp$cells)) && !all(temp$cells=="NA") && !all(temp$cells=="0")
      } else {flipflag <- FALSE}

      if (flipflag)
      {
# CASE THAT BOTH GMNC AND OPTIMUM LMXC ARE LESS THAN 1
# (BUT FIRST, CHECK THAT TLMXC ISN'T ALREADY AS LARGE AS IT CAN GET)
        if (min(1/tgmnc,length(trnl))!=max(1, 1/tlmxc))
        {
          tlmxc <-
            optimize( function(invlmxc) {
                ifelse(!all(is.na(SubDivStrat(rownames = tcnl,
                                        colnames = trnl,
                                        distspec = t(tdm),
                                        max.cpt = min(1/tgmnc, length(trnl)),
                                        min.cpt = invlmxc,
                                        tolerance = .5,
                                        omit.fraction = NULL)$cells)),
                       invlmxc, -invlmxc)
                },
                upper = min(1/tgmnc,length(trnl)),
                lower = max(1, 1/tlmxc),
                tol = 1,
                maximum = TRUE
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
                ifelse(!all(is.na(SubDivStrat( rownames=rown1, colnames=coln1, distspec=dist1,
                      min.cpt=max(gmnc1, 1/length(rown1)), max.cpt=lmxc1,
                      tolerance=.5, omit.fraction= switch(1+is.na(omf), omf,
                        NULL) )$cells)),
                  lmxc1, 2*length(coln1) - lmxc1)
                },
                lower=max(tgmnc,1), upper=min(length(tcnl), tlmxc), tol=1,
                maximum=FALSE, rown1=trnl, coln1=tcnl,
                dist1=tdm, gmnc1=tgmnc, omf1=omf
                )$objective )
        } }

    }
  return(tlmxc)
  }

  max.controls <- mapply(sps, min.controls, FUN = .maxControlsCap)

  return(list(given.min.controls = min.controls, strictest.feasible.max.controls = max.controls))
}
