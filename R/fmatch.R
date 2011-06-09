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

controls <- unique(dists$control)
treateds <- unique(dists$treatment)

if ((sum(dim(distances)[1]) + length(controls) + length(treateds)) > (1e+7-2))
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

nt <- length(treateds)
nc <- length(controls)
d <- dists$distance + 1
startn <- rep(c(1:nt), nc)
endn <- nt + rep(c(1:nc), rep(nt ,nc))
ucap <- rep(1, (nt*nc))

# Add figures for "end" and "sink" nodes
# "end" is node nt+nc+1; "sink" is node nt+nc+2
d <- c(d, rep(0, nc+nt), rep(0, nc))
startn <- c(startn, 1:(nt+nc), nt+ 1:nc)
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
return(x)

}
