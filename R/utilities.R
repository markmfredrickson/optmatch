################################################################################
# Utility functions used throughout the package
################################################################################

### toZ: turn various objects into a logical vector indicating treatment
### accepts a variety of inputs and keeps names/rownames if present

setGeneric("toZ", function(x) {
  if (is.data.frame(x) | is.matrix(x)) {
    if (dim(x)[2] > 1) {
      stop("Treatment indicators must be vectors or single column
      matrices/data.frames")
    }

    nms <- rownames(x)
    x <- x[,1]
    names(x) <- nms
  }

  nms <- names(x)
  tmp <- standardGeneric("toZ")
  names(tmp) <- nms
  return(tmp)
})

setMethod("toZ", "logical", function(x) {
  u <- unique(x)
  if (!(TRUE %in% u)) {
    stop("There must be at least one treatment unit.")
  }
  if (!(FALSE %in% u)) {
    stop("There must be at least one control unit.")
  }
  x
})

setMethod("toZ", "numeric", function(x) {
  u <- unique(x)
  if (any(!(u[!is.na(u)] %in% 0:1))) {
    stop("Numeric treatment indicators can only take on values 1 (treatment) and 0 (control).")
  }
  if (!(1 %in% u)) {
    stop("There must be at least one treatment unit.")
  }
  if (!(0 %in% u)) {
    stop("There must be at least one control unit.")
  }
  as.logical(x)
})

setOldClass("haven_labelled")
setMethod("toZ", "haven_labelled", function(x) {
  toZ(as.numeric(x))
})

setMethod("toZ", "character", function(x) {
  stop(paste("Character treatment indicators no longer supported.\n",
             "Convert into a numeric or logical vector."))
})

setMethod("toZ", "factor", function(x) {
  stop(paste("Factor treatment indicators no longer supported.\n",
             "Convert into a numeric or logical vector."))
})

#' (Internal) Remove the call before digesting a distance so things
#' like omitting caliper and calling caliper=NULL give the same digest
#' #
#' @param dist Distance object to hash. Must be one of
#'   \code{InfinitySparseMatrix}, \code{BlockedInfinitySparseMatrix},
#'   \code{DenseMatrix}, \code{matrix}, or \code{distmatch.dlist}.
#' @return Hash on the distance object with a null \code{call}
#' @keywords internal
#' @rdname dist_digest
dist_digest <- function(dist) {
  UseMethod("dist_digest")
}

#' @rdname dist_digest
dist_digest.DenseMatrix <- function(dist)  {
  dist_digest(dist@.Data)
}

#' @rdname dist_digest
dist_digest.matrix <- function(dist)  {
  sum(dist[is.finite(dist)]) + nrow(dist) + ncol(dist)
}

#' @rdname dist_digest
dist_digest.InfinitySparseMatrix <- function(dist) {
  sum(dist@.Data) + sum(dist@rows) + sum(dist@cols)
}

#' @rdname dist_digest
dist_digest.BlockedInfinitySparseMatrix <- function(dist) {
  dist_digest(as.InfinitySparseMatrix(dist))
}

#' @rdname dist_digest
dist_digest.optmatch.dlist <- function(dist) {
  # to do!
  1
}

#' (Internal) If the x argument does not exist for
#' match_on, fullmatch, or pairmatch, use this function
#' to print a helpful message.
#'
#' @param x_str x as a string of code, usually deparse(substitute(x))
#' @param data_str data arg as string of code
#' @param ... will look for 'z = <stuff>' in the extra args of caller
#' @return string a helpful error message
#' @author Josh Buckner
#' @keywords internal
missing_x_msg <- function(x_str, data_str, ...) {
  if(data_str == "NULL")
    data_str <- "<data argument>"

  # use a way of finding z that doesn't require z
  # to exist as an actual R obj
  z <- ""
  extra_args_str <- deparse(substitute(list(...)))
  z_regex <- "z = (\\S+)[,\\)]"
  z_search <- regexpr(z_regex, extra_args_str, perl=TRUE)
  if(z_search != -1) {
    z_match <- regmatches(extra_args_str, z_search)[1]
    z <- sub(z_regex, "\\1", z_match, perl=TRUE)
  }

  msg_tail <- if(z != "")
                paste("or ", z, "~", x_str, sep="")
              else
                ""

  paste("Can't find",
        paste(x_str, ".", sep=""),
        "If it lives within the data frame provided",
        "as the data argument, try",
        paste(data_str, "$", x_str, sep=""),
        msg_tail)
}
