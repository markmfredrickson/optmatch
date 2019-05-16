fmatch <- function(distance, max.row.units, max.col.units,
			min.col.units = 1, f = 1, stability.increment=1L, node_prices = NULL)
{
  if(!inherits(distance, "data.frame") & !all(colnames("data.frame") %in% c("treated", "control", "distance"))) {
    stop("Distance argument is not a canonical matching problem (an adjacency list of the graph): A data.frame with columns `treated`, `control`, `distance`.")
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
  treated.units <- levels(distance$treated)
  control.units <- levels(distance$control)
  nt <- length(treated.units)
  nc <- length(control.units)
  narcs <- dim(distance)[1]
  problem.size <- narcs + nt + nc

  if (problem.size > getMaxProblemSize()) {
      stop(paste('Maximum matching problem may have only',
                 signif(getMaxProblemSize(), 2), "- (nrows + ncols + 2) finite entries;",
                 problem.size - getMaxProblemSize(), 'too many.',
                 "Set 'options(\"optmatch_max_problem_size\" = Inf)' to disable this check."),
           call. = FALSE)
  }


  if (any(as.integer(distance$distance) != distance$distance |
	    distance$distance < 0)) {
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
  #warnings to prohibit use of names reserved for the two terminal nodes
  if (any(control.units == '(_Sink_)'))
  {
    warning('Cannot chose "(_Sink_)" or "(_End_)" as unit name')
  }
  if (any(control.units == '(_End_)'))
  {
    warning('Cannot chose "(_Sink_)" or "(_End_)" as unit name')
  }
  if (any(treated.units == '(_Sink_)'))
  {
    warning('Cannot chose "(_Sink_)" or "(_End_)" as unit name')
  }
  if (any(treated.units == '(_End_)'))
  {
    warning('Cannot chose "(_Sink_)" or "(_End_)" as unit name')
  }

  # set up the problem for the Fortran algorithm
  # each node has a integer ID number
  # startn indicates where each arc starts (using ID num)
  # endn indicates where each arc ends (using ID num)
  # nodes 1:nt are the treated units
  # nodes (nt + 1):nc are the control units
  # we use the levels of the treated and control factors to generate the ID numbers
  # the capacity of these arcs is 1

  dists <- as.vector(distance$distance) + stability.increment
  startn <- as.numeric(distance$treated)
  endn <- nt + as.numeric(distance$control)
  ucap <- rep(1, narcs)

  # Add entries for "end" and "sink" nodes
  # "end" is node nt+nc+1; "sink" is node nt+nc+2
  dists <- c(dists, rep(0, nc + nt), rep(0, nc))
  startn <- c(startn, 1:(nt + nc), nt + 1:nc)
  endn <- c(endn, rep(nt + nc + 1, nc + nt), rep(nt + nc + 2, nc))
  ucap <- c(ucap, rep(mxc - mnc, nt), rep(mxr - 1, nc), rep(1, nc))

  # supply
  b <- c(rep(mxc, nt), rep(0, nc), -(mxc * nt - round(f * nc)), -round(f * nc))

  if(!is.null(node_prices))
  {
    end.controls <- data.frame(control = "(_End_)", treated = treated.units, distance = 0)
    end.treatments <- data.frame(control = control.units, treated = "(_End_)", distance = 0)
    sink.control <- data.frame(control = control.units, treated = "(_Sink_)", distance = 0)
    distance <- rbind(distance, end.controls, end.treatments, sink.control)
    rcs <- prep.reduced.costs(distance, node_prices, narcs, nt, nc)
  }
  else
  {
    rcs <- as.integer(dists)
  }
  # If the user specifies using the old version of the relax algorithm. The `if` will be
  # FALSE if use_fallback_optmatch_solver is anything but TRUE, including NULL.
  # We have to duplicate the .Fortran code to make R CMD Check not complain about "registration" problems
  if (identical(options()$use_fallback_optmatch_solver, TRUE)) {
    fop <- .Fortran("relaxalgold",
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
  } else {
    fop <- .Fortran("relaxalg",
                    as.integer(nc + nt + 2),
                    as.integer(length(startn)),
                    as.integer(startn),
                    as.integer(endn),
                    as.integer(dists),
                    as.integer(ucap),
                    as.integer(b),
                    x1=integer(length(startn)),
                    rc1 = rcs,
                    crash1=as.integer(0),
                    large1=as.integer(.Machine$integer.max/4),
                    feasible1=integer(1),
                    NAOK = FALSE,
                    DUP = TRUE,
                    PACKAGE = "optmatch")
  }

  #browser()
  feas <- fop$feasible1 & ((mnc*nt <= round(f*nc) & mxc*nt >= round(f*nc)) |
            (round(f*nc) <= nt & round(f*nc)*mxr >= nt))

  x <- feas * fop$x1 - (1 - feas)

  if(is.null(node_prices))
  {
    end.controls <- data.frame(control = "(_End_)", treated = treated.units, distance = 0)
    end.treatments <- data.frame(control = control.units, treated = "(_End_)", distance = 0)
    sink.control <- data.frame(control = control.units, treated = "(_Sink_)", distance = 0)
    distance <- rbind(distance, end.controls, end.treatments, sink.control)
  }

  if (identical(options()$use_fallback_optmatch_solver, TRUE)) {
    ans <- x[1:narcs]
    rcosts <- fop$rc[1:narcs]
    cbind(distance[1:length(ans), ], solution = ans)
  } else
  {
    ans <- c(x[1:narcs], integer(length(fop$rc) - narcs))
    rcosts <- fop$rc
    obj <-cbind(distance, solution = ans, reduced.cost=rcosts)
    #sinkn.price <- fop$rc[which(startn == nt + 1 & endn == nt + nc + 2)] - fop$rc[which(startn == nt + 1 & endn == nt + nc + 1)]
    #attr(obj, 'sink.node.price') <- sinkn.price
    return(obj)
    #cbind(distance, solution = ans, reduced.cost=rcosts)

  }
}
#' @param df data.frame object containing all combinations of control:treated pairs and distances between them
#' @param node.prices Vector of node prices from a previously solved problem. At this point, these should be integer and adjusted to the current problems resolution
#' @param narcs.no.sink.or.end Number of arcs in the (sub)problem, excluding arcs connecting to sink or end nodes (?)
#' @param nt int, number of treated units
#' @param nc int, number of control units
#' @details This function converts node price and arc information (from a previously solved problem via warm start arguments) and converts it into reduced cost data to be passed along to the solver for warm start/initialization purposes. output is a vector of reduced costs, integer precision
prep.reduced.costs <- function(df, node.prices, narcs.no.sink.or.end, nt, nc)
{

  reduced.costs = numeric(nrow(df))
  # reduced.costs <- df$distance + node.prices[df$control] - node.prices[df$treated]
  reduced.costs[1:narcs.no.sink.or.end] <- df$distance[1:narcs.no.sink.or.end] + node.prices[as.character(df$control[1:narcs.no.sink.or.end])] - node.prices[as.character(df$treated[1:narcs.no.sink.or.end])]


  reduced.costs[(narcs.no.sink.or.end + 1):(nrow(df) - nc)] <- -c(node.prices[1:(nt + nc)])
  reduced.costs[(nrow(df)-nc +1):nrow(df)] <- node.prices["(_Sink_)"] - node.prices[(nt+1):(nt + nc)]
  return(as.integer(reduced.costs))
}
