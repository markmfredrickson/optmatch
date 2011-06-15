fmatch <- function(distances, max.row.units, max.col.units, 
			min.col.units=1, f=1)
{
###
# NB: ORDER OF ARGUMENTS SWITCHED FROM PREV VERSION!
mxc <- round(max.col.units) #  (formerly kt)
mnc <- round(min.col.units) #  (formerly ktl)
mxr <- round(max.row.units)
if (mnc>1) {mxr <- 1}

# Check that matching problem is well-specified
if (mxc<mnc) stop("min.col.units may not exceed max.col.units")
if (any(c(mxc, mnc, mxr)<1)) 
   { stop("max and min constraints must be 1 or greater") }
if (!is.numeric(f) | f>1 | f<0) stop("f must be a number in [0,1]")

stopifnot(is(distances, "DistanceSpecification"))
dists <- prepareMatching(distances)
# dists is a data.frame with 3 columns: control, treatment, distance
# for defensive reasons, we make sure the control and treatment columns are factors
dists$control <- as.factor(dists$control) ; dists$treatment <- as.factor(dists$treatment)

# the Fortan code uses numeric ids. the first nt + nc are resevered for the units
# to be matched. 1:nt are treated (nt + 1):(nt + nc) are the controls
treateds <- as.integer(dists$treatment)
nt <- length(unique(treateds))

controls <- as.integer(dists$control) + nt 
nc <- length(unique(controls))

if ((sum(dim(distances)[1]) + nc + nt) > (1e+7-2))
  stop(paste('matrix arg to fmatch may have only',
             1e+7, "-(nrows+ncols+2) finite entries;",
             sum(finiteind) + sum(dim(distance.matrix)) - 1e+7 -2, 'too many'),
       call. = FALSE)

if (any(any(dists$distances < 0)))
  stop("Distances must be positive")

# S4 TODO: Promote non-integer values into integers instead of an error
# if (any(as.integer(distance.matrix[finiteind])!=distance.matrix[finiteind] | 
# 	distance.matrix[finiteind]<0)) 
#    { stop("distance.matrix should be nonnegative integer") }
# SHOULD PROBABLY DISABLE NEXT TWO WARNINGS
if (mxc!=max.col.units | mnc!=min.col.units | 
	(mxr==round(max.row.units) & mxr!=max.row.units) ) 
   { warning("fmatch coerced one or more constraints to integer") }
if (mnc>1 & round(max.row.units)>1) 
   { warning("since min.col.units>1, fmatch coerced max.row.units to 1") }

d <- dists$distance
startn <- treateds 
endn <- controls
ucap <- rep(1, (nt*nc))

# Add figures for "end" and "sink" nodes
# "end" is node nt+nc+1; "sink" is node nt+nc+2
d <- c(d, rep(0, nc+nt), rep(0, nc))
startn <- c(startn, 1:(nt+nc), nt+ 1:nc) # NOTE: this assumes treates and controls form a 1:(nc+nt) vector
endn <- c(endn, rep(nt+nc+1, nc+nt), rep(nt+nc+2, nc))
ucap <- c(ucap, rep(mxc-mnc, nt), rep(mxr-1, nc), rep(1, nc))

# supply
b <- c(rep(mxc, nt), rep(0, nc), -(mxc*nt-round(f*nc)), -round(f*nc))

fop <- .Fortran("relaxalg", 
		as.integer(nc+nt+2), 
		as.integer(length(startn)), 
		as.integer(startn), 
		as.integer(endn), 
		as.integer(d), 
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

x <- feas*fop$x - (1-feas)
# trim the pseudo nodes (e.g. the sink) from the list
x <- x[1:(dim(dists)[1])]

return(matches2factor(dists, x))

}

# a helper for fmatch(),
# turns the cannonical form and the results of the Fortran code into a factor of matches
# NOTE: this algorithm _should_ run in linear time. It rather depends on how good data.frame look up
# is. If that is itself linear (and not constant), then this algorithm runs in quadratic time, which would be
# nice to avoid.
matches2factor <- function(data, matches) {
  matches.only <- data[matches == 1, ]
  controls <- levels(as.factor(data$control))
  treateds <- levels(as.factor(data$treatment))

  # for each treated units, get the controls it is linked to and vice-versa
  treated.control <- lapply(treateds, function(t) { matches.only$control[matches.only$treatment == t]})
  names(treated.control) <- treateds
  control.treated <- lapply(controls, function(t) { matches.only$treatment[matches.only$control == t]})
  names(control.treated) <- controls

  # now we see which treateds are linked to the same group of controls.
  # the groups vector will be a set of numerics, each indicating a matched group.
  nt <- length(treated.control)
  groupst <- vector("numeric", nt)
  for (i in nt:1) { # counts down to give better numbering in the sets.
    val <- treated.control[[i]]
    same <- as.logical(lapply(treated.control, function(i) { identical(i, val) }))
    groupst[same] <- i
  }
  names(groupst) <- treateds
  
  # now we look up the control ids from the treated ids
  nc <- length(control.treated)
  groupsc <- vector("numeric", nt)
  for (i in nc:1) { # counts down to give better numbering in the sets.
    tied.to <- control.treated[[i]][1]
    groupsc[i] <- groupst[tied.to]
  }
  names(groupsc) <- controls

  tmp <- as.factor(paste("m", c(groupst, groupsc), sep = ""))
  names(tmp) <- c(treateds, controls)
  return(tmp)  
}
