fmatch <- function(distance, max.row.units, max.col.units, 
			min.col.units = 1, f = 1)
{
  # distance must have a prepareMatching object
  if (!hasMethod("prepareMatching", class(distance))) {
    stop("Argument \'distance\' must have a \'prepareMatching\' method")  
  }
  # NB: ORDER OF ARGUMENTS SWITCHED FROM PREV VERSION!
  mxc <- round(max.col.units) #  (formerly kt)
  mnc <- round(min.col.units) #  (formerly ktl)
  mxr <- round(max.row.units)
  
  if (mnc > 1) {
    mxr <- 1
  }

  # Check that matching problem is well-specified
  if (mxc < mnc) {
    stop("min.col.units may not exceed max.col.units")
  } 

  if (any(c(mxc, mnc, mxr) < 1)) { 
    stop("max and min constraints must be 1 or greater") 
  }

  if (!is.numeric(f) | f > 1 | f < 0) {
    stop("f must be a number in [0,1]")
  } 

  # Ok, let's see how a big a problem we are working on
  # this returns a "canonical" matching: a data.frame with
  # three columns: treated, control, distance. Where the first two
  # are factors and the last is numeric.
  distance.prepared <- prepareMatching(distance)
  treated.units <- levels(distance.prepared$treated)
  control.units <- levels(distance.prepared$control)
  nt <- length(treated.units)
  nc <- length(control.units) 
  narcs <- dim(distance.prepared)[1]
  problem.size <- narcs + nt + nc

  # these "soft" limits are backed by "hard" limits in the Fortran code itself
  if (problem.size > (1e+7-2)) {
      stop(paste('matrix arg to fmatch may have only',
                 1e+7, "-(nrows+ncols+2) finite entries;",
                 problem.size - 1e+7 - 2, 'too many'),
           call. = FALSE)
  }

  
  if (any(as.integer(distance.prepared$distance) != distance.prepared$distance | 
	    distance.prepared$distance < 0)) { 
    stop("distance should be nonnegative integer") 
  }

  # SHOULD PROBABLY DISABLE NEXT TWO WARNINGS
  if (mxc != max.col.units | mnc!=min.col.units | 
  	  (mxr == round(max.row.units) & mxr != max.row.units)) { 
      warning("fmatch coerced one or more constraints to integer")
  }

  if (mnc > 1 & round(max.row.units) > 1) { 
    warning("since min.col.units > 1, fmatch coerced max.row.units to 1") 
  }

  # set up the problem for the Fortan algorithm
  # each node has a integer ID number
  # startn indicates where each arc starts (using ID num)
  # endn indicates where each arc ends (using ID num)
  # nodes 1:nt are the treated units
  # nodes (nt + 1):nc are the control units
  # we use the levels of the treated and control factors to generate the ID numbers
  # the capacity of these arcs is 1

  dists <- as.vector(distance.prepared$distance) + 1
  startn <- as.numeric(distance.prepared$treated)
  endn <- nt + as.numeric(distance.prepared$control)
  ucap <- rep(1, narcs)

  # Add entries for "end" and "sink" nodes
  # "end" is node nt+nc+1; "sink" is node nt+nc+2
  dists <- c(dists, rep(0, nc + nt), rep(0, nc))
  startn <- c(startn, 1:(nt + nc), nt + 1:nc)
  endn <- c(endn, rep(nt + nc + 1, nc + nt), rep(nt + nc + 2, nc))
  ucap <- c(ucap, rep(mxc - mnc, nt), rep(mxr - 1, nc), rep(1, nc))

  # supply
  b <- c(rep(mxc, nt), rep(0, nc), -(mxc * nt - round(f * nc)), -round(f * nc))

  fop <- .Fortran("relaxalg", 
		as.integer(nc + nt + 2), 
		as.integer(length(startn)), 
		as.integer(startn), 
		as.integer(endn), 
		as.integer(dists), 
		as.integer(ucap), 
		as.integer(b), 
		x1=integer(length(startn)), 
		crash1=as.integer(0), 
		large1=as.integer(.Machine$integer.max/4), 
		feasible1=integer(1),
		NAOK = FALSE,
		DUP = TRUE,
		PACKAGE = "optmatch")

  feas <- fop$feasible & ((mnc*nt <= round(f*nc) & mxc*nt >= round(f*nc)) | 
            (round(f*nc) <= nt & round(f*nc)*mxr >= nt))

  x <- feas * fop$x - (1 - feas)

  ans <- numeric(narcs)
  ans <- x[1:narcs]
  return(cbind(distance.prepared, solution = ans))

  ### Not evaluated becuase we are passing more information back to the caller.
  # this using the helper below, 
  res <- numeric(sum(dim(distance)))
  names(res) <- c(rownames(distance), colnames(distance))
  res[names(tmp)] <- tmp

  return(res)
}


