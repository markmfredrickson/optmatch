################################################################################
# Utility functions used throughout the package
################################################################################

### toZ: turn various objects into a logical vector indicating treatment
### accepts a variety of inputs and keeps names/rownames if present

setGeneric("toZ", function(x) {

  if (any(is.na(x))) {
    stop("NAs not allowed in treatment indicator.")
  }

  if (is.data.frame(x) | is.matrix(x)) {
    if (dim(x)[2] > 1) {
      stop("Treatment indicators must be vectors or single column
      matrices/data.frames")
    }

    nms <- rownames(x)
    x <- x[,1]
    names(x) <- nms
  }

  if (length(unique(x)) != 2) {
    stop(paste("Treatment indicator must have exactly 2 levels not",
      length(unique(x))))
  }

  nms <- names(x)
  tmp <- standardGeneric("toZ")
  names(tmp) <- nms
  return(tmp)
})

# a noop
setMethod("toZ", "logical", identity)

# we already know it has two levels, so just call as.logical
setMethod("toZ", "numeric", function(x) as.logical(x))

setMethod("toZ", "character", function(x) toZ(as.factor(x)))

setMethod("toZ", "factor", function(x) {
  toZ(as.numeric(x) - 1)
})

#' (Internal) Remove the call before digesting a distance so things
#' like omitting caliper and calling caliper=NULL give the same digest
#' #
#' @param dist Distance object to hash. Must be one of
#'   \code{InfinitySparseMatrix}, \code{BlockedInfinitySparseMatrix},
#'   \code{DenseMatrix}, \code{matrix}, or \code{distmatch.dlist}.
#' @return Hash on the distance object with a null \code{call}
#' @import digest
dist_digest <- function(dist) {
  if (class(dist)[1] %in% c("InfinitySparseMatrix", "BlockedInfinitySparseMatrix", "optmatch.dlist", "DenseMatrix", "matrix")) {
    csave <- attr(dist, "call")
    attr(dist, "call") <- NULL
    out <- digest::digest(dist)
    attr(dist, "call") <- csave
    return(out)
  }
  stop("Must pass distance object")
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
