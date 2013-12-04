SubDivStrat <- function(rownames, colnames, distspec, min.cpt,
    max.cpt, tolerance, omit.fraction=NULL, matched.distances=FALSE)
{
  if (min.cpt <=0 | max.cpt<=0) {
    stop("inputs min.cpt, max.cpt must be positive")
  }

  if (!all(rownames %in% dimnames(distspec)[[1]])) {
    stop("input \'rownames\' may only contain row names of input \'distspec\'")
  }

  if (!all(colnames %in% dimnames(distspec)[[2]])) {
    stop("input \'rownames\' may only contain col. names of input \'distspec\'")
  }

  # distance must have a prepareMatching object
  if (!hasMethod("prepareMatching", class(distspec))) {
    stop("Argument \'distspec\' must have a \'prepareMatching\' method")
  }

  # convert the distspec into a cannonical matching specification with columns
  # treated, control, distance
  dm <- prepareMatching(distspec)

  rownames <- as.character(rownames)
  colnames <- as.character(colnames)

  rfeas <- length(unique(dm$treated))
  cfeas <- length(unique(dm$control))

  # If any controls were unmatchable, they were dropped by prepareMatching, and omit.fraction
  # needs to be updated.
  if (cfeas < length(colnames) & is.numeric(omit.fraction)) {
    omit.fraction <- (omit.fraction*length(colnames) - (length(colnames) - cfeas))/cfeas
    # This can happen if the number to be omitted is less than the number of unmatchables
    if (omit.fraction <= 0) {
      omit.fraction <- NULL
    }
  }

  if (floor(min.cpt) > ceiling(max.cpt) | ceiling(1/min.cpt) < floor(1/max.cpt))
  {
    ans <- rep("NA",length(rownames)+length(colnames))
      names(ans) <- c(rownames, colnames)
      return(list(cells=ans, maxerr=NULL, distance=NULL))
  }

  # the next block of code, the dm <- ... is commented out as
  # dm is no longer a matrix. Completely unreachable entries may be a
  # problem later, but
  if (is.null(omit.fraction)) {
    f.ctls <- 1
    # dm <- matrix(dm[rfeas, cfeas], sum(rfeas), sum(cfeas),
    # dimnames=list(rownames[rfeas], colnames[cfeas]))
  } else {
    if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction > 1) {
      stop("omit.fraction must be null or between 0 and 1")
    }

    f.ctls <- 1-omit.fraction
    # dm <- matrix(dm[rfeas,], sum(rfeas), length(colnames),
    #              dimnames=list(rownames[rfeas], colnames))
  }

  if (any(rfeas) & any(cfeas))
  {
    old.o <- options(warn=-1)

    if (any(dm$distance > 0)) {
      reso <- (.Machine$integer.max/64 -2)/max(dm$distance)
    } else {
      reso <- min(.Machine$integer.max/64 -2, (sum(rfeas)+sum(cfeas))/tolerance)
    }

    if (tolerance>0 & sum(rfeas)>1 & sum(cfeas)>1) {
      reso <- min(reso, (sum(rfeas) + sum(cfeas) - 2)/tolerance)
    }

    options(old.o)

    # fmatch returns a matrix with columns `treatment`, `control`, and `solution`
    # it also has a column `distance` with toIntFuction(dm * reso)
    .matcher <- function(toIntFunction) {
      tmp <- dm
      tmp$distance <- toIntFunction(dm$distance * reso)
      fmatch(tmp,
             max.row.units = ceiling(1/min.cpt),
             max.col.units = ceiling(max.cpt),
             min.col.units = max(1, floor(min.cpt)), f=f.ctls)

    }

    temp <- .matcher(floor)

    if (any(is.na(temp$solution))) {
      maxerr <- 0
    } else {
      maxerr <- sum(temp$solution * dm$distance, na.rm = TRUE) -
                sum(temp$solution * temp$distance, na.rm = TRUE) / reso +
                (sum(rfeas) > 1 & sum(cfeas) > 1) *
                (sum(rfeas) + sum(cfeas) - 2 - sum(temp$solution)) / reso
    }

    if (maxerr > tolerance)
    {
      temp1 <- temp
      temp2 <- .matcher(round)

      if  (sum(temp1$solution * dm$distance, na.rm = TRUE) <= sum(temp2$solution * dm$distance, na.rm = TRUE)) {
        temp <- temp1
      } else {
        temp <- temp2
      }

      maxerr <- sum(temp$solution * dm$distance, na.rm = TRUE) -
                sum(temp1$solution * temp$distance, na.rm = TRUE)/reso +
                (max(1, sum(rfeas) - 1) + max(1, sum(cfeas) - 1) -
                  (sum(rfeas) == 1 & sum(cfeas) == 1) - sum(temp1$solution)) / reso
    }

    ### NOTE: this if statment not yet updated to new fmatch data format.
    if (matched.distances) {
      if (all(!is.na(temp$solution))) {
        dma <- max(dm[as.logical(temp$solution)])
        dist <- c(apply(temp$solution * pmin(dm, dma), 1, sum),
                  apply(temp$solution * pmin(dm, dma), 2, sum))

        dist[c(rep(FALSE, dim(temp)[1]),
               apply(temp * apply(temp, 1, sum), 2, sum) == 1)] <- NA

        dist[c(apply(temp, 1, sum) > 1, apply(temp, 2, sum) > 1)] <- NA

      } else{
        dist <- rep(NA, sum(dim(temp)))
        mode(dist) <- "numeric"
        names(dist) <- unlist(dimnames(temp))
      }
    } else {
      dist <- 0
    }
  } else {
    temp <- 0 ; maxerr <- 0 ; dist <- 0
  }

  ans <- rep(NA,length(rownames)+length(colnames))
  names(ans) <- c(rownames, colnames)

  matches <- solution2factor(temp)
  ans[names(matches)] <- matches

  return(list(cells = ans, err = maxerr))
}

# a small helper function to turn a solution data.frame into a factor of matches
solution2factor <- function(s) {
  s2 <- s[s$solution == 1,]

  if (dim(s2)[1] == 0) {
    return(NULL)
  }

  # control units are labeled by the first treated unit to which they are connected
  # unlist(as.list(...)) was the best way I could find to make this into a vector, keeping names
  control.links <- unlist(as.list(by(s2, s2$control, function(x) { x[1,"treated"] })))

  # treated units are labeld by the label of the first control unit to which they are connected
  treated.links <- unlist(as.list(by(s2, s2$treated, function(x) { control.links[x[1, "control"]][1] })))

  # join the links
  return(c(treated.links, control.links))

}
