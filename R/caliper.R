#' @include InfinitySparseMatrix.R
#' @include DenseMatrix.R
NA

#' Prepare matching distances suitable for matching within calipers.
#'
#' Encodes calipers, or maximum allowable distances within which to
#' match. The result of a call to \code{caliper} is itself a distance specification between
#' treated and control units that can be used with 
#' \code{pairmatch()} or \code{fullmatch()}. Calipers can also be combined with
#' other distance specifications for richer matching problems.
#'
#' \code{caliper} is a generic function with methods for any of the allowed distance
#' specifications: user created matrices, the results of \code{\link{match_on}}, the results
#' of \code{\link{exactMatch}}, or combinations (using \code{`+`}) of these objects.
#'
#' \code{width} provides the size of the caliper, the allowable distance for
#' matching. If the distance between a treated and control pair is less than or equal
#' to this distance, it is allowed kept; otherwise, the pair is discarded from future
#' matching. The default comparison of ``equal or less than can'' be changed to any other
#' comparison function using the \code{comparison} argument.
#' 
#' If you wish to exclude specific units from the caliper requirements, pass the names of 
#' these units in the \code{exclude} argument. These units will be allowed to match any other
#' unit.
#'
#' @author Mark M. Fredrickson and Ben B. Hansen
#'
#' @references  P.~R. Rosenbaum and D.~B. Rubin (1985),
#' \sQuote{Constructing a control group using multivariate matched sampling
#'   methods that incorporate the propensity score},
#'  \emph{The American Statistician}, \bold{39} 33--38.
#'
#' @seealso \code{\link{exactMatch}}, \code{\link{match_on}}, \code{\link{fullmatch}}, \code{\link{pairmatch}}
#'
#' @example inst/examples/caliper.R
#'
#' @keywords nonparametric
#' 
#' @param x A distance specification created with \code{\link{match_on}} or similar.
#' @param width The width of the caliper: how wide of a margin to
#' allow in matches. Be careful in setting the width. Vector valued
#' arguments will be recycled for each of the finite entries in
#' \code{x} (and no order is guaranteed for \code{x} for some types of
#' distance objects).
#' @param exclude (Optional) A character vector of observations (corresponding to row and column names) to exclude from the caliper.
#' @param compare A function that decides that whether two
#' observations are with the caliper. The default is
#' \code{`<=`}. \code{`<`} is a common alternative.
#' @param values Should the returned object be made of all zeros
#' (\code{values = FALSE}, the default) or should the object include
#' the values of the original object (\code{values = TRUE})?
#' @return A matrix like object that is suitable to be given
#' as \code{distance} argument to \code{\link{fullmatch}} or
#' \code{\link{pairmatch}}. The caliper will be only zeros and \code{Inf} values,
#' indicating a possible match or no possible match, respectively.
#'
#' You can combine the results of \code{caliper} with other distances using the
#' \code{`+`} operator. See the examples for usage.
#' @export
#' @docType methods
#' @rdname caliper-methods
setGeneric("caliper", function(x, width = 1, exclude = c(), compare = `<=`, values = FALSE) {
  tmp <- standardGeneric("caliper")
  tmp@call <- match.call()
  return(tmp)
})

#' @rdname caliper-methods
#' @aliases caliper,InfinitySparseMatrix-method
setMethod("caliper", "InfinitySparseMatrix",
function(x, width = 1, exclude = c(), compare = `<=`, values = FALSE) {

  excluded.rows <- which(x@rownames %in% exclude)
  excluded.cols <- which(x@colnames %in% exclude)

  y <- discardOthers(x, as.vector(compare(x, width)) | 
                     x@rows %in% excluded.rows |
                     x@cols %in% excluded.cols)

  if (!values) {
    y@.Data <- rep(0, length(y@.Data))
  }

  return(y)
})

#' @rdname caliper-methods
#' @aliases caliper,matrix-method
setMethod("caliper", "matrix",
function(x, width = 1, exclude = c(), compare = `<=`, values = FALSE) {
  caliper(as.InfinitySparseMatrix(x), width = width, exclude = exclude, compare = compare, values = values)  
})

#' @rdname caliper-methods
#' @aliases caliper,optmatch.dlist-method
setMethod("caliper", "optmatch.dlist",
function(x, width = 1, exclude = c(), compare = `<=`, values = FALSE) {
  caliper(as.matrix(x), width = width, exclude = exclude, compare = compare, values = values)  
})
