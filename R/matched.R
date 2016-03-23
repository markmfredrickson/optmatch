##' Given a bipartite matching (object of class \code{optmatch}),
##' create a logical vector of the same length indicating which units
##' were and were not placed into matched sets.
##'
##' \code{matched} and \code{unmatched} indicate which elements of
##' \code{x} do and do not belong to matched sets, as indicated by
##' their character representations in \code{x}.
##'
##' When \code{fullmatch} has been presented with an inconsistent
##' combination of constraints and discrepancies between potential
##' matches, so that there exists no matching (i) with finite total
##' discrepancy within matched sets that (ii) respects the given
##' constraints, then the matching problem is said to be infeasible.
##' \code{TRUE}s in the output of \code{matchfailed} indicate that
##' this has occurred.
##'
##' @note To understand the output of \code{matchfailed} element-wise,
##'   note that \code{fullmatch} handles a matching problem in three
##'   steps.  First, if \code{fullmatch} has been directed to match
##'   within subclasses, then it divides its matching problem into a
##'   subproblem for each subclass.  Second, \code{fullmatch} removes
##'   from each subproblem those individual units that lack
##'   permissible potential matches (i.e. potential matches from which
##'   they are separated by a finite discrepancy).  Such ``isolated''
##'   units are flagged in such a way as to be indicated by
##'   \code{unmatched}, but not by \code{matchfailed}.  Third,
##'   \code{fullmatch} presents each subproblem, with isolated
##'   elements removed, to an optimal matching routine.  If such a
##'   reduced subproblem is found at this stage to be infeasible, then
##'   each unit contributing to it is so flagged as to be indicated by
##'   \code{matchfailed}.
##'
##' @title Identification of units placed into matched sets
##'
##' @param x Vector of class \code{optmatch} (especially as generated
##'   by a call to \code{fullmatch}).
##'
##' @return A logical vector (without names).
##'
##' @seealso \code{\link{fullmatch}}
##'
##' @examples
##' data(plantdist)
##'
##' mxpl.fm0 <- fullmatch(plantdist) # A feasible matching problem
##' c(sum(matched(mxpl.fm0)), sum(unmatched(mxpl.fm0)))
##' sum(matchfailed(mxpl.fm0))
##' mxpl.fm1 <- fullmatch(plantdist, # An infeasible problem
##'                       max.controls=3, min.controls=3)
##' c(sum(matched(mxpl.fm1)), sum(unmatched(mxpl.fm1)))
##' sum(matchfailed(mxpl.fm1))
##'
##' mxpl.si <- factor(c('a', 'a', 'c', rep('d',4), 'b', 'c', 'c', rep('d', 16)))
##' names(mxpl.si) <- LETTERS[1:26]
##' mxpl.exactmatch <- exactMatch(mxpl.si, c(rep(1, 7), rep(0, 26 - 7)))
##' # Subclass a contains two treated units but no controls;
##' # subclass b contains only a control unit;
##' # subclass c contains one treated and two control units;
##' # subclass d contains the remaining twenty units.
##' # only valid subproblems will be used
##'
##' mcl <- c(1, Inf)
##'
##' mxpl.fm2 <- fullmatch(plantdist + mxpl.exactmatch,
##'                       max.controls=mcl)
##' sum(matched(mxpl.fm2))
##'
##' table(unmatched(mxpl.fm2), matchfailed(mxpl.fm2))
##'
##' mxpl.fm2[matchfailed(mxpl.fm2)]
##'
##' mxpl.fm2[unmatched(mxpl.fm2) &   # isolated units return as
##'          !matchfailed(mxpl.fm2)] # unmatched but not matchfailed
##'
##' @keywords manip
##' @author Ben Hansen
##' @export
##' @rdname matched
matched <- function(x) !is.na(x)

#' @rdname matched
#' @export
unmatched <- is.na
