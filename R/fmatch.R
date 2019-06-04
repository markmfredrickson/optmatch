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
##* @param node_info NodeInfo object for this subproblem, or `NULL` 
##* @return data frame with columns `treatment`, `control`, `solution`, `distance`. 
##*         The `distance` column holds values of the integer distance that `intSolve()` saw. 
##* @author Ben Hansen, Mark Fredrickson, Adam Rauh
##* @keywords internal

fmatch <- function(distance, max.row.units, max.col.units,
			min.col.units = 1, f = 1, node_info = NULL)
{
    if (identical(options()$use_fallback_optmatch_solver, TRUE))
       warning("Old version of RELAX-IV solver (avoiding variable-sized Fortran arrays)\n no longer implemented; using current version.")
    if (!inherits(distance, "data.frame") ||
        !setequal(colnames(distance), c("control", "treated", "distance"))
        ) 
        stop("Distance argument is not a canonical matching problem\n (an adjacency list of the graph): A data.frame with columns `treated`, `control`, `distance`.")
  
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
      if (any(distance[['distance']] > .Machine$integer.max)) stop("Integer distances too large for machine.\n (Tolerance/resolution too fine?)")
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

  ## distance is a "canonical" matching: a data.frame with
  ## three columns, 'control', 'treated' and 'distance', the first two
  ## being factors and the last being a numeric.  Typically
  ## 'control' and 'treated' are interpretable as control (Z=0)
  ## and treated (Z=1) in the matching problem as it was presented
  ## to fullmatch(), pairmatch() or match_on(), but in full
  ## matching the subproblem may have been "flipped" prior  
  ## to its delivery to this function.  In that case, what's  
  ## "treated" here would be control in the originating problem,  
  ## and vice versa.  To eliminate ambiguities in such instances, 
  ## we'll prefer "row.units" or "upstream" to "treated",
  ## and also "col.units" or "downstream" over "control". 
  row.units <- levels(distance$treated)
  col.units <- levels(distance$control)
  nt <- length(row.units)
  nc <- length(col.units)
  n.mc  <- as.integer(round(nc * f)) # no. downstream (usu., control) units to be matched
  narcs <- nrow(distance)

  ## Ok, let's see how a big a problem we are working on    
  problem.size <- narcs + nt + 2L*nc # no. arcs in MCF rep'n of problem
  if (problem.size > getMaxProblemSize()) {
      stop(paste('Maximum matching problem may have only',
                 signif(getMaxProblemSize(), 2), "- (nrows + 2*ncols) finite entries;",
                 problem.size + 1L - getMaxProblemSize(), 'too many.',
                 "Set 'options(\"optmatch_max_problem_size\" = Inf)' to disable this check."),
           call. = FALSE)
  }

  
  if (mnc > 1 & mxr > 1) {
    stop("since min.col.units > 1, max.row.units can be at most 1.")
  }
  #prohibit use of names reserved for the two terminal nodes
  if (any(c(row.units, col.units) == '(_Sink_)'))
    stop('Cannot choose "(_Sink_)" as unit name.')
  if (any(c(row.units, col.units) == '(_End_)'))
    stop('Cannot choose "(_End_)" as unit name')

  ## Bypass solver if problem is recognizably infeasible
  if ( (mxr >1 & nt/mxr > n.mc) | #max.row.units too low
       (mxr==1L & nt * mnc > n.mc) |# min.col.units too high  
       (nt * mxc < n.mc) #max.col.units too low
      )
  {
    return(cbind(distance[1:narcs, ], solution = rep(-1L, narcs)))
      }

  ##  Min-Cost-Flow representation of problem  ####
  ## Each node has a integer ID number, implicitly pointing
  ## to corresponding row of this data frame.
  ## nodes 1:nt are the upstream matchable units in the flow diagram
  ## nodes (nt + 1):nc are the downstream matchable units, etc.
  nodes  <- new("NodeInfo",
            data.frame(name=c(row.units, col.units,
                              "(_End_)", "(_Sink_)"),# note that `factor()` would put these at start not end of levels
                       price=0L, 
                       upstream_not_down=c(rep(TRUE, nt), rep(FALSE, nc), NA, NA), 
                       supply=c(rep(mxc, nt), rep(0L, nc),
                                -(mxc * nt - n.mc), -n.mc),
                       groups=factor(rep(NA_character_, nt+nc+2L)),
                       stringsAsFactors=FALSE
                       )
            )
    if (!is.null(node_info))
        {
            stopifnot(is(node_info, "NodeInfo"),
                      is.integer(node_info$price),
                      !any(is.na(node_info$price)),
                      all(row.units %in% node_info[['name']]),
                      sum(node_info[['name']]=="(_End_)")==1,
                      sum(node_info[['name']]=="(_Sink_)")==1)
            nodes[,'price']  <- node_info[match(c(row.units, col.units,
                                                  "(_End_)", "(_Sink_)"),
                                                node_info$name),'price']
            ## if a col unit doesn't have a match w/in node_info$name, then
            ## if nodes now has an NA for its price. Replace that w/ min of
            ## bookkeeping node prices, so that CS has a chance of being preserved.
            if (length(col.noprice  <- setdiff(col.units, node_info[['name']])))
                nodes[nodes[['name']] %in% col.noprice, 'price']  <-
                    min(node_info[match(c("(_End_)", "(_Sink_)"), node_info$name),
                                  'price', drop=TRUE]
                        )
            }
  # ID numbers implicitly point to corresp. row of nodes
  ####    Arcs involving End or Sink nodes   ####
  bookkeeping  <-
      data.frame(groups=factor(rep(NA_character_, nt +2L*nc)),
                 start=c(1L:(nt + nc), # ~ c(row.units, col.units)
                         nt + 1L:nc), # ~ col.units
                 end=c(rep(nt + nc + 1L, nc + nt), # nt + nc+ 1L ~ '(_End_)'
                       rep(nt + nc + 2L, nc) ), # nt + nc + 2L ~ '(_Sink_)'   
          flow=0L,
          capacity=c(rep(mxc - mnc, nt), rep(mxr - 1L, nc), rep(1L, nc))
                 )
  ## set up the problem for the Fortran algorithm ##
  ## startn indicates where each arc starts (using ID num)
  ## endn indicates where each arc ends (using ID num)
  ## "End"/"Overflow" is node nt+nc+1; "Sink" is node nt+nc+2
  #################### All arcs  ###################
  ###         Arcs rep'ing potential matches    Arcs involving End or Sink
  startn<- c(as.integer(distance$treated),     bookkeeping$start )
  endn <-  c(nt + as.integer(distance$control),bookkeeping$end )
  ucap <-  c(rep(1L, narcs),                   bookkeeping$capacity )
  redcost<-c(distance$distance,                rep(0L, nrow(bookkeeping)) ) +
        nodes[endn, "price"] - nodes[startn, "price"]
    
  ## Everything headed for the solver ought currently to be cast as integer,
  ## but I couldn't immediately find documentation of rules for preserving integer
  ## cast in arithmetic ops, and it's hard to be sure those won't change.  Since
  ## feeding a double by mistake will cause a difficult to diagnose infeasibility,
  ## let's just be sure:
  if (!is.integer(problem.size)) problem.size  <- as.integer(problem.size)    
  if (!is.integer(startn)) startn  <- as.integer(startn)    
  if (!is.integer(endn)) endn  <- as.integer(endn)
  if (!is.integer(ucap)) ucap  <- as.integer(ucap)
  if (!is.integer(nodes$supply)) nodes$supply  <- as.integer(nodes$supply)
  if (!is.integer(redcost)) redcost  <- as.integer(redcost)  
  
    fop <- .Fortran("relaxalg",
                    n1=as.integer(nc + nt + 2L),
                    na1=problem.size,
                    startn1=startn,
                    endn1=endn,
                    c1=c(distance$distance, rep(0L, nrow(bookkeeping))),
                    u1=ucap,
                    b1=nodes$supply,
                    x1=integer(problem.size),
                    rc1 = redcost,
                    crash1=as.integer(0),
                    large1=as.integer(.Machine$integer.max/4),
                    feasible1=integer(1),
                    NAOK = FALSE,
                    DUP = TRUE,
                    PACKAGE = "optmatch")

  ## Material used to create s3 optmatch object:
  feas <- fop$feasible1
  x <- feas * fop$x1 - (1 - feas)
  obj <-cbind(distance, solution = x[1:narcs])

  #### Recover node prices, store in nodes table ##
  ## In full matching, each upstream (row) or downstream (column) node starts
  ## an arc ending at End, and these are also the only arcs 
  ## ending there. End being at the bottom of the canonical diagram, 
  ## call these "downarcs". Downarcs' costs are 0, so their reduced
  ## costs are simply the price of the End node minus the price  
  ## of the upstream or downstream node they started from. Imposing 
  ## a convention that the price of End is 0, the prices of these upstream
  ## and downstream nodes are just the opposites of the reduced
  ## costs of the arcs they begin.   
  End_rownum_in_nodes_table  <- which(nodes$name=="(_End_)")
  stopifnot(length(End_rownum_in_nodes_table)==1)
  is_downarc  <-  ( bookkeeping[['end']] == End_rownum_in_nodes_table )
  downarc_redcosts  <- fop$rc1[c(rep(FALSE, narcs), is_downarc )]  
  nodes[ bookkeeping[['start']][ is_downarc ] ,
                                      "price" ]  <-
      -1L * downarc_redcosts 
  ## It remains to recover the price of the Sink node.  There 
  ## are arcs to it from every downstream (column) node, and these
  ## are the only arcs involving it.  First we extract these arcs'  
  ## reduced costs.
  Sink_rownum_in_nodes_table  <- which(nodes$name=="(_Sink_)")
  stopifnot(length(Sink_rownum_in_nodes_table)==1)
  is_arctosink  <-  ( bookkeeping[['end']] == Sink_rownum_in_nodes_table )
  arctosink_redcosts  <- fop$rc1[c(rep(FALSE, narcs), is_arctosink )]
  arctosink_price_diffs  <-   
      nodes[ bookkeeping[['start']][ is_arctosink ] ,
                                            "price" ]
  sinkprice  <- arctosink_redcosts + arctosink_price_diffs
    if (any(is.na(sinkprice))) {
        sinkprice  <- as.double(arctosink_redcosts) + as.double(arctosink_price_diffs)
        if (any(abs(sinkprice - sinkprice[1])>.Machine$double.eps^.5))
            stop("Inferred sink prices not mutually consistent.")
            } else {
                if (!all(sinkprice==sinkprice[1]))
                    stop("Mutually inconsistent inferred sink prices.")
                }
  nodes[nodes$name=="(_Sink_)", "price"]  <- sinkprice[1]  

  ### Recover arc flow info, store in `arcs` ###
  ## info extracted from problem solution:
  matches  <- distance[as.logical(fop$x1[1L:narcs]), c("treated", "control")]
  bookkeeping[1L:(problem.size-narcs), "flow"]  <- fop$x1[(narcs+1L):problem.size]
  ## reshape `matches`, `bookkeeping` to match ArcInfo object spec:
  colnames(matches)  <- c("upstream", "downstream")
  matches  <- data.frame(groups=factor(rep(NA_integer_, nrow(matches))),
                         matches)  
  bookkeeping[['start']]  <- factor(nodes$name[ bookkeeping[['start']] ])
  bookkeeping[['end']]  <- factor(nodes$name[ bookkeeping[['end']] ])
  arcs  <- new("ArcInfo", matches=matches, bookkeeping=bookkeeping)

  sp  <- new("SubProbInfo")
  sp[1L, "feasible"]  <- TRUE
  fmcfs  <- new("FullmatchMCFSolutions", subproblems=sp,
                nodes=nodes, arcs=arcs, matchables=new("MatchablesInfo"))
    c(obj,
      list(maxerr=0), # if we were called from doubleSolve(), this will be re-set there
      list(MCFSolution=fmcfs) )
}
