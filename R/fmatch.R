##* This is the function that calls the solver
##*
##* Also handles pair matching, 1:k matching and matching
##* with between k>1 and l>k controls.  (If you want matching
##* with between k>1 and l>k members of the treatment group
##* per matched set, then you have to transpose the (virtual)
##* distance matrix before feeding the problem to this function.)
##*
##* Distances should be integer and **positive** (no 0s).  There
##* can be NAs: these will be converted to the minimum of provided
##* distances, unless all provided distances are NA, in which case
##* they'll be interpreted as 1L's. 
##* @title Full matching via RELAX-IV min cost flow solver
##* @param distance data frame w/ integer columns distance, treated, control; see Details
##* @param max.row.units numeric, upper limit on num treated units per matched set
##* @param max.col.units numeric, upper limit on num control units per matched set
##* @param min.col.units numeric, lower limit on num control units per matched set
##* @param f double, fraction of all row units to be matched
##* @param node_prices 
##* @return data frame with columns `treatment`, `control`, `solution`, `distance`. 
##*         The `distance` column holds values of the integer distance that `intSolve()` saw. 
##* @author Ben Hansen, Mark Fredrickson, Adam Rauh
##* @keywords internal

fmatch <- function(distance, max.row.units, max.col.units,
			min.col.units = 1, f = 1, node_prices = NULL)
{
  if(!inherits(distance, "data.frame") & !all(colnames("data.frame") %in% c("treated", "control", "distance"))) {
    stop("Distance argument is not a canonical matching problem (an adjacency list of the graph): A data.frame with columns `treated`, `control`, `distance`.")
  }
  stopifnot(is.numeric(f))
  mxc <- as.integer(round(max.col.units))
  mnc <- as.integer(round(min.col.units))
  mxr <- as.integer(round(max.row.units))

  if (mnc > 1) {
    mxr <- 1L
  }

  # Check that matching problem is well-specified
  if (mxc < mnc) {
    stop("min.col.units may not exceed max.col.units")
  }

  if (any(c(mxc, mnc, mxr) < 1)) {
    stop("max and min constraints must be 1 or greater")
  }

  if (!is.integer(distance[['distance']])) {
      tdist  <- as.integer(distance[['distance']])
      if (isFALSE(all.equal(distance[['distance']],tdist))) stop("distance should be integer")
      distance[['distance']]  <- tdist
  }
  if (any(nadists  <-  is.na(distance[['distance']])))
  {
      replacement  <- if (all(nadists)) { #occurs if
                          1L              #user passed ISM 
                      } else {            #produced by exactMatch()`
                          min(distance[['distance']], na.rm=TRUE)
                          }
      distance[['distance']][nadists]  <- replacement
  }
  
  if (mxr > 1) # i.e. many-one matches permissible
  {
      if (any(distance[['distance']] <= 0))
          stop("Nonpositive discrepancies not compatible with full matching\n (see Hansen & Klopfer, 2006 JCGS, sec.4.1).")
      } else if (any(distance[['distance']] < 0)) stop("distance should be nonnegative")

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
  narcs <- nrow(distance)
  problem.size <- narcs + nt + nc
  n.mc  <- as.integer(round(nc * f)) # number of matched controls

  if (problem.size > getMaxProblemSize()) {
      stop(paste('Maximum matching problem may have only',
                 signif(getMaxProblemSize(), 2), "- (nrows + ncols + 2) finite entries;",
                 problem.size - getMaxProblemSize(), 'too many.',
                 "Set 'options(\"optmatch_max_problem_size\" = Inf)' to disable this check."),
           call. = FALSE)
  }

  
  if (mnc > 1 & mxr > 1) {
    stop("since min.col.units > 1, max.row.units can be at most 1.")
  }
  #prohibit use of names reserved for the two terminal nodes
  if (any(c(treated.units, control.units) == '(_Sink_)'))
    stop('Cannot choose "(_Sink_)" as unit name.')
  if (any(c(treated.units, control.units) == '(_End_)'))
    stop('Cannot choose "(_End_)" as unit name')

  ## Bypass solver if problem is recognizably infeasible
  if ( (mxr >1 & nt/mxr > n.mc) | #max.row.units too low
       (mxr==1L & nt * mnc > n.mc) |# min.col.units too high  
       (nt * mxc < n.mc) #max.col.units too low
      )
  {
    return(cbind(distance[1:narcs, ], solution = rep(-1L, narcs)))
      }
  
  # set up the problem for the Fortran algorithm
  # each node has a integer ID number
  # startn indicates where each arc starts (using ID num)
  # endn indicates where each arc ends (using ID num)
  # nodes 1:nt are the treated units
  # nodes (nt + 1):nc are the control units
  # we use the levels of the treated and control factors to generate the ID numbers
  # the capacity of these arcs is 1

  dists <- as.vector(distance$distance) 
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
  b <- c(rep(mxc, nt), rep(0, nc), -(mxc * nt - n.mc), -n.mc)

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

  feas <- fop$feasible1

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
    return(obj)
  }
}

prep.reduced.costs <- function(df, node.prices, narcs.no.sink.or.end, nt, nc)
{
  reduced.costs = numeric(nrow(df))
  # reduced.costs <- df$distance + node.prices[df$control] - node.prices[df$treated]
  reduced.costs[1:narcs.no.sink.or.end] <- df$distance[1:narcs.no.sink.or.end] + node.prices[as.character(df$control[1:narcs.no.sink.or.end])] - node.prices[as.character(df$treated[1:narcs.no.sink.or.end])]


  reduced.costs[(narcs.no.sink.or.end + 1):(nrow(df) - nc)] <- -c(node.prices[1:(nt + nc)])
  reduced.costs[(nrow(df)-nc +1):nrow(df)] <- node.prices["(_Sink_)"] - node.prices[(nt+1):(nt + nc)]
  return(as.integer(reduced.costs))
}
