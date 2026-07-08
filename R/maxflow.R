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
##* satisfiable if and only if the max flow value is
##* min.controls * (number of row units); row units whose source arcs
##* remain unsaturated in that computation are reported as deficient.
##* With source arcs capped at max.controls, the max flow value is the
##* largest number of column units usable in any matching respecting the
##* per-set maximum.  Design assumption (to be verified against the
##* solver in tests): when the lower bounds are satisfiable, that maximum
##* is attainable jointly with them, so every total in
##* [min.controls * nt, V] is feasible.
##*
##* @title Max-flow feasibility check for matching without control sharing
##* @param distspec distance specification with an `edgelist` method
##*        (matrix, InfinitySparseMatrix, ...); rows are treated, columns control
##* @param min.controls numeric of length 1, at least 1; minimum number of
##*        column units per row unit
##* @param max.controls numeric of length 1, at least min.controls; maximum
##*        number of column units per row unit
##* @return list with elements `feasible` (can every row unit
##*         simultaneously receive min.controls column units?),
##*         `max.controls.matchable` (max flow value: largest number of
##*         column units usable), `min.omit.fraction` and
##*         `max.mean.controls` (translations of that value; NA if not
##*         feasible), and `deficient.row.units` (character; row units
##*         left short of min.controls in the particular max flow found —
##*         max flows are not unique, so this attribution is one witness,
##*         not canonical)
##* @author Josh Buckner, Ben Hansen
##* @keywords internal
maxflow_feasibility <- function(distspec, min.controls = 1,
                                max.controls = min.controls) {
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

  lower <- flow.value(mnc)
  feasible <- (lower[["cost"]] == mnc * nt)
  deficient <- row.units[lower[["flows"]][1L:nt] < mnc]

  v <- if (mxc == mnc) lower[["cost"]] else flow.value(mxc)[["cost"]]

  list(feasible = feasible,
       max.controls.matchable = as.integer(v),
       min.omit.fraction = if (feasible) (nc - v) / nc else NA_real_,
       max.mean.controls = if (feasible) v / nt else NA_real_,
       deficient.row.units = deficient)
}
