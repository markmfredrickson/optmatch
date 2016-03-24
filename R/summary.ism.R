##' Summarize a distance matrix
##'
##' Given a distance matrix, return information above it, including
##' dimension, sparsity information, unmatchable members, summary of
##' finite distances, and, in the case of
##' \code{BlockedInfinitySparseMatrix}, block structure.
##'
##' The output consists of several pieces.
##'
##' \itemize{
##'  \item{Membership: }{Indicates the dimension of the distance.}
##'
##'  \item{Total (in)eligible potential matches: }{A measure of the
##'   sparsity of the distance. Eligible matches have a finite
##'   distance between treatment and control members; they could be
##'   matched. Ineligible matches have \code{Inf} distance and can not
##'   be matched. A higher number of ineligible matches can speed up
##'   matching, but runs the risk of less optimal overall matching
##'   results.}
##'
##'  \item{Unmatchable treatment/control members: }{If any
##'   observations have no eligible matches (e.g. their distance to
##'   every potential match is \code{Inf}) they are listed here. See
##'   Value below for details of how to access lists of matchable and
##'   unmatchable treatment and control members.}
##'
##'  \item{Summary of minimum matchable distance per treatment member:
##'   }{To assist with choosing a caliper, this is a numeric summary
##'   of the smallest distance per matchable treatment member. If you
##'   provide a caliper that is less than the maximum value, at least
##'   one treatment member will become unmatchable.}
##'
##'  \item{Block structure: }{For \code{BlockedInfinitySparseMatrix},
##'   a quick summary of the structure of each individual block. (The
##'   above will all be across all blocks.) This may indicate which
##'   blocks, if any, are problematic.}
##' }
##'
##' @param object A \code{InfinitySparseMatrix},
##'   \code{BlockedInfinitySparseMatrix} or \code{DenseMatrix}.
##' @param ... Ignored.
##' @param distanceSummary Default \code{TRUE}. Should a summary of
##'   minimum distance per treatment member be calculated? May be slow
##'   on larger data sets.
##' @param printAllBlocks If \code{object} is a
##'   \code{BlockedInfinitySparseMatrix}, should summaries of all
##'   blocks be printed alongside the overall summary? Default
##'   \code{FALSE}.
##' @param blockStructure If \code{object} is a
##'   \code{BlockedInfinitySparseMatrix} and \code{printAllBlocks} is
##'   false, print a quick summary of each individual block. Default
##'   \code{TRUE}. If the number of blocks is high, consider
##'   suppressing this.
##' @return A named \code{list}. The summary for an
##'   \code{InfinitySparseMatrix} or \code{DenseMatrix} contains the
##'   following:
##'
##'   \itemize{
##'
##'    \item{\code{total}: }{Contains the total number of treatment
##'     and control members, as well as eligible and ineligible
##'     matches.}
##'
##'    \item{\code{matchable}: }{The names of all treatment and
##'     control members with at least one eligible match.}
##'
##'    \item{\code{unmatchable}: }{The names of all treatment and
##'     control members with no eligible matches.}
##'
##'    \item{\code{distances}: }{The summary of minimum matchable
##'     distances, if \code{distanceSummary} is \code{TRUE}.}
##'
##'   }
##'
##'   For \code{BlockedInfinitySparseMatrix}, the named \code{list}
##'   instead of contains one entry per block, named after each block
##'   (i.e. the value of the blocking variable) as well as a block
##'   named 'overall' which contains the summary ignoring blocks. Each
##'   of these entries contains a \code{list} with entries 'total',
##'   'matchable', 'unmatchable' and 'distances', as described above.
##' @export
##' @name summary.ism
summary.InfinitySparseMatrix <- function(object, ..., distanceSummary=TRUE) {

  finitedata <- is.finite(object@.Data)
  rowsfinite <- object@rows[finitedata]
  colsfinite <- object@cols[finitedata]
  datafinite <- object@.Data[finitedata]

  mtreat <- 1:dim(object)[1] %in% sort(unique(rowsfinite))
  mcontrol  <- 1:dim(object)[2] %in% sort(unique(colsfinite))

  if (distanceSummary & length(datafinite)) {
    distances <- summary(sapply(split(datafinite, rowsfinite), min))
  } else {
    distances <- NULL
  }

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)
  attr(out, "ismname") <- deparse(substitute(object))

  class(out) <- "summary.InfinitySparseMatrix"
  out
}

##' @export
##' @rdname summary.ism
summary.BlockedInfinitySparseMatrix <- function(object, ...,
                                                distanceSummary=TRUE,
                                                printAllBlocks=FALSE,
                                                blockStructure=TRUE) {

  ismname <- deparse(substitute(object))

  out <- lapply(levels(object@groups),
                function(x) {
                  thisgroup <- names(object@groups[object@groups == x])
                  ism <- subset(object,
                                subset=object@rownames %in% thisgroup,
                                select=object@colnames %in% thisgroup)
                  s <- summary(ism, ..., distanceSummary=distanceSummary)
                  attr(s, "ismname") <- ismname
                  attr(s, "blockname") <- x
                  return(s)
                })
  names(out) <- levels(object@groups)

  out$overall <- summary.InfinitySparseMatrix(object, ...,
                                              distanceSummary=distanceSummary)

  attr(out, "ismname") <- ismname
  attr(out, "blocknames") <- levels(object@groups)

  attr(out$overall, "ismname") <- attr(out, "ismname")

  attr(out, "printAllBlocks") <- printAllBlocks
  attr(out, "blockStructure") <- blockStructure

  class(out) <- "summary.BlockedInfinitySparseMatrix"
  return(out)
}

##' @export
##' @rdname summary.ism
summary.DenseMatrix <- function(object, ..., distanceSummary=TRUE) {
  mtreat <- apply(object, 1, function(x) any(is.finite(x)))
  mcontrol <- apply(object, 2, function(x) any(is.finite(x)))
  if (distanceSummary & length(object@.Data[is.finite(object@.Data)])) {
    distances <- summary(apply(object, 1, min))
  } else {
    distances <- NULL
  }

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)

  attr(out, "ismname") <- deparse(substitute(object))

  class(out) <- "summary.DenseMatrix"
  out
}

internal.summary.helper <- function(x,
                                    matchabletxt,
                                    matchablectl,
                                    distances=NULL) {
  out <- list()
  d <- dim(x)

  # Size of treatment and control groups
  out$total$treatment <- d[1]
  out$total$control <- d[2]

  # Count of eligble and ineligible pairs.
  out$total$matchable <- Reduce("+", num_eligible_matches(x))
  out$total$unmatchable <- prod(d) - out$total$matchable

  out$matchable$treatment <- rownames(x)[matchabletxt]
  out$matchable$control <- colnames(x)[matchablectl]
  out$unmatchable$treatment <- rownames(x)[!matchabletxt]
  out$unmatchable$control <- colnames(x)[!matchablectl]

  out$distances <- distances
  return(out)
}

##' @export
print.summary.InfinitySparseMatrix <- function(x, ...) {
  cat(paste("Membership:", x$total$treatment, "treatment,",
            x$total$control, "control\n"))
  cat(paste("Total eligible potential matches:", x$total$matchable,
            "\n"))
  cat(paste("Total ineligible potential matches:", x$total$unmatchable,
            "\n"))
  cat("\n")

  numunmatch <- sapply(x$unmatchable, length)
  for (i in 1:2) {
    if (numunmatch[i] > 0) {
      cat(paste0(numunmatch[i], " unmatchable ", names(numunmatch)[i],
                 " member", if(numunmatch[i] > 1) { "s" } , ":\n"))
      cat("\t")
      cat(paste(x$unmatchable[[i]][1:min(5, numunmatch[i])],
                collapse=", "))
      if (numunmatch[i] > 5) {
        cat(", ...\n")

        cat(paste0("See summary(", attr(x, "ismname"), ")",
                   if (!is.null(attr(x, "blockname"))) {
                     paste0("$`", attr(x, "blockname"), "`")
                   }, "$unmatchable$",
                   names(numunmatch)[i], " for a complete list."))
        }
      cat("\n\n")
    }
  }

  if (!is.null(x$distances) && any(!is.na(x$distances))) {
    cat("Summary of minimum matchable distance per treatment member:\n")
    print(x$distances, ...)
    cat("\n")
  }
  return(invisible(x))
}

##' @export
print.summary.BlockedInfinitySparseMatrix <- function(x, ...) {

  cat("Summary across all blocks:\n")
  print(x$overall, ...)
  blockentries <- names(x) %in% attr(x, "blocknames")

  if (!attr(x, "printAllBlocks")) {
    if (attr(x, "blockStructure")) {
      cat("Block structure:\n")
      blocksummary <- matrix(unlist(lapply(x[blockentries], "[", "total")),
                             byrow=TRUE, ncol=4)
      blocksummary <- cbind(sapply(sapply(sapply(x[blockentries], "[", "matchable"), "[", "treatment"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "matchable"), "[", "control"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "unmatchable"), "[", "treatment"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "unmatchable"), "[", "control"), length))
      rownames(blocksummary) <- paste0("`",attr(x, "blocknames"),"`")
      colnames(blocksummary) <- c("Matchable Txt",
                                  "Matchable Ctl",
                                  "Unmatchable Txt",
                                  "Unmatchable Ctl")
      print(blocksummary)

      cat("\n")
    }

    cat(paste0("To see summaries for individual blocks,",
               " call for example summary(",
               attr(x, "ismname"), ")$`",
               attr(x, "blocknames")[1], "`.\n"))
  } else {
    cat("Indiviual blocks:\n\n")
    for (i in attr(x, "blocknames")) {
      cat(paste0("`",i,"`\n"))
      print(x[[i]])
    }
  }

  cat("\n")
  return(invisible(x))
}

##' @export
print.summary.DenseMatrix <- function(x, ...) {
  print.summary.InfinitySparseMatrix(x, ...)
  return(invisible(x))
}
