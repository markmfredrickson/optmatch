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
##* @param distance EdgeList; see details
##* @param max.row.units numeric, upper limit on num treated units per matched set
##* @param max.col.units numeric, upper limit on num control units per matched set
##* @param min.col.units numeric, lower limit on num control units per matched set
##* @param f double, fraction of all row units to be matched
##* @param node_info NodeInfo object for this subproblem
##* @return data frame with columns `dist`, `i`, `j`, `solution`.
##*         The `dist` column holds values of the integer distance that `intSolve()` saw.
##* @author Ben Hansen, Mark Fredrickson, Adam Rauh
##* @keywords internal

fmatch <- function(distance, max.row.units, max.col.units,
			min.col.units = 1, f = 1, node_info)
{
    if (identical(options()$use_fallback_optmatch_solver, TRUE))
       warning("Old version of RELAX-IV solver (avoiding variable-sized Fortran arrays)\n no longer implemented; using current version.")
    stopifnot(is(distance, "EdgeList"), is(node_info, "NodeInfo"),
              is.numeric(f))
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

  thenodes  <- node_info[['name']]
  row.units <- subset(thenodes, node_info[['upstream_not_down']])
  col.units <- subset(thenodes, !node_info[['upstream_not_down']])
    distance  <-
        edgelist(distance,
                 c(row.units,# NB: solver prep below assumes row.units come
                   col.units)# before col.units in levels(distance$i/j)
                 )# This also removes arcs to/from nodes not in subproblem.
  if (!is.integer(distance[['dist']])) {
      tdist  <- as.integer(distance[['dist']])
      if (isFALSE(all.equal(distance[['dist']],tdist))) stop("distance should be integer")
      distance  <- edgelist(data.frame(i=distance[['i']], j=distance[['j']],
                                       dist=tdist),
                            c(row.units, col.units)
                            )
  }
  if (any(nadists  <-  is.na(distance[['dist']])))
  {
      replacement  <- if (all(nadists)) { #occurs if
                          1L              #user passed ISM 
                      } else {            #produced by exactMatch()`
                          min(distance[['dist']], na.rm=TRUE)
                          }
      distance[['dist']][nadists]  <- replacement
  }
  
  if (mxr > 1) # i.e. many-one matches permissible
  {
      if (any(distance[['dist']] <= 0))
          stop("Nonpositive discrepancies not compatible with full matching\n (see Hansen & Klopfer, 2006 JCGS, sec.4.1).")
      } else if (any(distance[['dist']] < 0)) stop("distance should be nonnegative")

  if (!is.numeric(f) | f > 1 | f < 0) {
    stop("f must be a number in [0,1]")
  }

  ## distance is an EdgeList rep'n of a "canonical" matching, with
  ## columns, 'j', 'i' and 'dist', the first two
  ## being factors and the last being a numeric.  Typically
  ## 'j' and 'i' are interpretable as control (Z=0)
  ## and treated (Z=1) in the matching problem as it was presented
  ## to fullmatch(), pairmatch() or match_on(), but in full
  ## matching the subproblem may have been "flipped" prior  
  ## to its delivery to this function.  In that case, what's  
  ## "treated" here would be control in the originating problem,  
  ## and vice versa.  To eliminate ambiguities in such instances, 
  ## we'll prefer "i', "row.units" or "upstream" to "treated",
  ## and also "j", "col.units" or "downstream" over "control".
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
    node.labels(nodes)  <- nodes[['name']] # new convention per i166

    if ( !all(node_info[['price']]==0) )
    {
        nodes[,'price']  <- node_info[match(c(row.units, col.units,
                                              "(_End_)", "(_Sink_)"),
                                            node_info$name),
                                      'price']
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
  startn<- c(as.integer(distance[['i']]), bookkeeping$start )
  endn <-  c(as.integer(distance[['j']]), bookkeeping$end )
  ucap <-  c(rep(1L, narcs),              bookkeeping$capacity )
  redcost<-c(distance[['dist']],               rep(0L, nrow(bookkeeping)) ) +
        nodes[endn, "price"] - nodes[startn, "price"]
    
  ## Everything headed for the solver ought currently to be cast as integer,
  ## but passing a double by mistake will cause a difficult-to-diagnose
  ## infeasibility.  So let's just be sure:
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
                    c1=c(distance[['dist']], rep(0L, nrow(bookkeeping))),
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
  sinkprice  <- arctosink_redcosts +
      nodes[ bookkeeping[['start']][ is_arctosink ] ,
                                            "price" ]
  if (!all(sinkprice==sinkprice[1])) stop("Mutually inconsistent inferred sink prices.")
  nodes[nodes$name=="(_Sink_)", "price"]  <- sinkprice[1]  

  ### Recover arc flow info, store in `arcs` ###
  ## info extracted from problem solution:
  matches  <- distance[as.logical(fop$x1[1L:narcs]), c("i", "j")]
  bookkeeping[1L:(problem.size-narcs), "flow"]  <- fop$x1[(narcs+1L):problem.size]
  ## reshape `matches`, `bookkeeping` to match ArcInfo object spec:
  matches  <- data.frame(groups=factor(rep(NA_integer_, nrow(matches))),
                         upstream=factor(matches[['i']],
                                         levels=node.labels(nodes)
                                         ),
                         downstream=factor(matches[['j']],
                                         levels=node.labels(nodes)
                                         )
                         )
  bookkeeping[['start']]  <- factor(bookkeeping[['start']], 
                                    levels=1L:(nt + nc + 2),
                                    labels=node.labels(nodes)
                                    )
  bookkeeping[['end']]  <- factor(bookkeeping[['end']],
                                  levels=1L:(nt + nc + 2),
                                  labels=node.labels(nodes)
                                  )
  arcs  <- new("ArcInfo", matches=matches, bookkeeping=bookkeeping)

  sp  <- new("SubProbInfo")
  sp[1L, "feasible"]  <- TRUE
  fmcfs  <- new("FullmatchMCFSolutions", subproblems=sp,
                nodes=nodes, arcs=arcs)
    c(obj,
      list(maxerr=0), # if we were called from doubleSolve(), this will be re-set there
      list(MCFSolution=fmcfs) )
}
