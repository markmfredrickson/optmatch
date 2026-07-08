##* Max-flow feasibility assessment for matching without control sharing (#200)
##*
##* Prototype for issue #200: instead of settling for some feasible
##* omit.fraction during infeasibility recovery, compute the smallest
##* feasible omit.fraction (equivalently, the largest feasible
##* mean.controls) by solving a max flow problem via rlemon.
##*
##* This handles only problems without sharing of controls, i.e.
##* min.controls >= 1: pair matching, 1:k matching and matching with
##* between k and l>=k controls per treated unit.  The network is the
##* bipartite one: a source node with an arc of capacity max.controls
##* (resp. min.controls, see below) to each row unit, an arc of capacity 1
##* from row unit i to column unit j wherever distspec records a finite
##* discrepancy, and an arc of capacity 1 from each column unit to a sink
##* node.  Distances play no role beyond marking which pairings are
##* permitted, so no discretization is needed.
##*
##* Two max flow computations are made.  With source arcs capped at
##* min.controls, the problem's lower bounds are simultaneously
##* satisfiable if and only if the max flow value is min.controls times
##* the number of row units required to be matched; row units whose
##* source arcs remain unsaturated in that computation are reported as
##* deficient.  With source arcs capped at max.controls, the max flow
##* value is the largest number of column units usable in any matching
##* respecting the per-set maximum.  Design assumption (verified against
##* the solver in tests): when the lower bounds are satisfiable, that
##* maximum is attainable jointly with them, so every total in
##* [min.controls * (rows required), max flow value] is feasible.
##*
##* Row units without any finite discrepancy can never be matched.  With
##* drop.isolated.rows = FALSE they render the problem infeasible; with
##* drop.isolated.rows = TRUE they are exempted from the min.controls
##* requirement and reported in `unmatchable.row.units`, mirroring
##* solve_reg_fm_prob()'s practice of dropping such units from the node
##* table and leaving them unmatched.
##*
##* @title Max-flow feasibility check for matching without control sharing
##* @param distspec distance specification with an `edgelist` method
##*        (matrix, InfinitySparseMatrix, ...); rows are treated, columns control
##* @param min.controls numeric of length 1, at least 1; minimum number of
##*        column units per row unit
##* @param max.controls numeric of length 1, at least min.controls; maximum
##*        number of column units per row unit
##* @param drop.isolated.rows logical; exempt row units with no finite
##*        discrepancies from the min.controls requirement?
##* @return list with elements `feasible` (can every required row unit
##*         simultaneously receive min.controls column units?),
##*         `max.controls.matchable` (max flow value: largest number of
##*         column units usable), `min.omit.fraction` and
##*         `max.mean.controls` (translations of that value; NA if not
##*         feasible), `deficient.row.units` (character; required row
##*         units left short of min.controls in the particular max flow
##*         found -- max flows are not unique, so this attribution is one
##*         witness, not canonical), and `unmatchable.row.units`
##*         (character; row units with no finite discrepancies at all)
##* @author Josh Buckner, Ben Hansen
##* @keywords internal
maxflow_feasibility <- function(distspec, min.controls = 1,
                                max.controls = min.controls,
                                drop.isolated.rows = FALSE) {
  mnc <- as.integer(round(min.controls))
  mxc <- as.integer(round(max.controls))
  if (mnc < 1)
    stop("min.controls must be at least 1 (no sharing of controls)")
  if (mxc < mnc)
    stop("min.controls may not exceed max.controls")

  row.units <- dimnames(distspec)[[1]]
  col.units <- dimnames(distspec)[[2]]
  nt <- length(row.units)
  nc <- length(col.units)
  stopifnot(nt > 0, nc > 0)

  dm <- edgelist(distspec, c(row.units, col.units))
  ## Node IDs: 1:nt row units, nt + (1:nc) column units (the ordering
  ## edgelist() gives factor codes), then source and sink.
  keep <- is.finite(dm[["dist"]])
  match.arc.i <- as.integer(dm[["i"]])[keep]
  match.arc.j <- as.integer(dm[["j"]])[keep]
  narcs <- length(match.arc.i)

  isolated <- setdiff(1L:nt, match.arc.i)
  required <- if (drop.isolated.rows) setdiff(1L:nt, isolated) else 1L:nt
  n.req <- length(required)

  sourceID <- nt + nc + 1L
  sinkID <- nt + nc + 2L
  arc.sources <- c(rep(sourceID, nt), match.arc.i, nt + 1L:nc)
  arc.targets <- c(1L:nt, match.arc.j, rep(sinkID, nc))

  flow.value <- function(cap.per.row) {
    capacities <- c(rep(cap.per.row, nt), rep(1L, narcs), rep(1L, nc))
    rlemon::MaxFlow(arcSources = arc.sources,
                    arcTargets = arc.targets,
                    arcCapacities = capacities,
                    sourceNode = sourceID,
                    destNode = sinkID,
                    numNodes = sinkID)
  }

  ## isolated rows carry no flow, so requiring cost == mnc * n.req is
  ## exactly the simultaneous-lower-bounds condition on the required rows
  lower <- flow.value(mnc)
  feasible <- (lower[["cost"]] == mnc * n.req)
  deficient <- row.units[required[lower[["flows"]][required] < mnc]]

  v <- if (mxc == mnc) lower[["cost"]] else flow.value(mxc)[["cost"]]

  list(feasible = feasible,
       max.controls.matchable = as.integer(v),
       min.omit.fraction = if (feasible) (nc - v) / nc else NA_real_,
       max.mean.controls = if (feasible && n.req > 0) v / n.req else NA_real_,
       deficient.row.units = deficient,
       unmatchable.row.units = row.units[isolated])
}

##* Compose a human-readable explanation of why omitting controls cannot
##* make a no-sharing subproblem feasible, from a maxflow_feasibility()
##* result.  Used for #226-style messaging in fullmatch()'s recovery.
##* @param mf list as returned by maxflow_feasibility()
##* @param min.controls the subproblem's minimum controls per treated unit
##* @param n.asked number of controls the user's constraints ask to match
##* @return character of length 1
##* @keywords internal
maxflow_infeasibility_message <- function(mf, min.controls, n.asked) {
  listsome <- function(units) {
    if (length(units) > 5)
      paste0(paste(units[1:5], collapse = ", "), ", ... (",
             length(units), " units in all)")
    else
      paste(units, collapse = ", ")
  }
  msgs <- character(0)
  if (length(mf$unmatchable.row.units) > 0)
    msgs <- c(msgs, paste0("treated unit(s) with no permissible control: ",
                           listsome(mf$unmatchable.row.units)))
  if (!mf$feasible)
    msgs <- c(msgs, paste0("treated unit(s) that cannot be given ",
                           "min.controls (=", min.controls, ") controls ",
                           "without taking controls from other treated ",
                           "units, e.g.: ",
                           listsome(mf$deficient.row.units)))
  else if (mf$max.controls.matchable == 0)
    msgs <- c(msgs, "no treated unit has a permissible control")
  else if (n.asked <= mf$max.controls.matchable)
    msgs <- c(msgs, paste0("the requested omit.fraction leaves only ",
                           n.asked, " controls to be matched, too few to ",
                           "give every matchable treated unit min.controls ",
                           "(=", min.controls, ") controls"))
  paste0("A subproblem is infeasible and omitting further controls cannot ",
         "make it feasible: ", paste(msgs, collapse = "; "), ".")
}
