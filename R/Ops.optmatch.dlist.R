#' @export
Ops.optmatch.dlist <- function (e1, e2=NULL)
{
    ok <- switch(.Generic, "%%" = , "%/%" = FALSE, TRUE)
    if (!ok) {
        warning(.Generic, " not meaningful for matching distances; returning 1st arg")
        return(e1)
    }

    unary <- nargs() == 1

    if (nchar(.Method[1])) {
     rn1 <- attr(e1, "row.names")
     nne <- unlist(as.logical(lapply(e1, length)))
     e1.nullentries <- e1[!nne]
     full.sc1 <- names(e1)
     e1 <- e1[nne]
     sc1 <- names(e1)
     } else {rn1 <- NULL}

   if (nchar(.Method[2])) {
     rn2 <- attr(e2, "row.names")
     nne <- unlist(as.logical(lapply(e2, length)))
     e2.nullentries <- e2[!nne]
     full.sc2 <- names(e2)
     e2 <- e2[nne]
     sc2 <- names(e2)
     } else {rn2 <- NULL}

    if (!unary && all(nchar(.Method)))
      {
    rn12rn2 <- match(rn1, rn2)
    rn22rn1 <- match(rn2, rn1)
    if (any(is.na(rn12rn2)) && any(is.na(rn22rn1))) stop("arguments\' row names attributes don't match")
    if (!any(is.na(rn12rn2)) && any(diff(rn12rn2)<0)) stop("arguments\' row names inconsistently ordered")
    if (!any(is.na(rn22rn1)) && any(diff(rn22rn1)<0)) stop("arguments\' row names inconsistently ordered")

    # the proper behavior is:
    # - make sure the two objects have same length
    # - in each item, make sure the row and column names are the same
    # if either is not met, fail

    if (setequal(sc1,sc2)) {
      # if they have the same names, great. proceed, perhaps reording e2
      e2 <- e2[sc1]
    } else {
      if (length(sc1) != length(sc2)) {
        stop("arguments must have equal number of subproblems")
      }

      k <- length(sc1)
      for (i in 1:k) {
        if (!identical(dimnames(e1[[i]]), dimnames(e2[[i]]))) {
          stop("arguments must have identically named subproblem matrices")
        }
      }
    }


     dm11 <- lapply(e1, function(x) {if (is.null(dim(x))) {length(x)} else {dim(x)[1]}})
     dm11 <- unlist(dm11)
     dm12 <- lapply(e2, function(x) {if (is.null(dim(x))) {1} else {dim(x)[2]}})
     dm12 <- unlist(dm12)
     dm21 <- lapply(e2, function(x) {if (is.null(dim(x))) {length(x)} else {dim(x)[1]}})
     dm21 <- unlist(dm21)
     dm22 <- lapply(e2, function(x) {if (is.null(dim(x))) {1} else {dim(x)[2]}})
     dm22 <- unlist(dm22)

     if (any(dm11!=dm21) || any(dm12!=dm22))
       stop("dimensions of distlist arguments don\'t match")
  }
    value <- list()
    FUN <- get(.Generic, envir = parent.frame(), mode = "function")
     f <- if (unary)
        quote(FUN(left))
    else quote(FUN(left, right))



    if (nchar(.Method[1]) )
      {
      for (j in 1:length(e1))
        {
          left <- e1[[j]]
          if (!unary) {
            if (nchar(.Method[2])) {
              right <- e2[[j]] } else {
                right <- e2}
        }
          value[[j]] <- eval(f)
        }

      names(value) <- sc1

      if (length(e1.nullentries))
        {
          value <- c(value, e1.nullentries)
          value <- value[full.sc1]
        }
    } else
    {
      if (nchar(.Method[2]))
        {
      for (j in 1:length(e2))
        {
          right <- e2[[j]]
          left <- e1
          value[[j]] <- eval(f)
        }

      names(value) <- sc2

          if (length(e2.nullentries))
            {
            value <- c(value, e2.nullentries)
            value <- value[full.sc2]
            }
        }
    }

    class(value) <- c('optmatch.dlist', 'list')
    if (length(rn1)>length(rn2))
      {
        attr(value, "row.names") <- rn1
        } else {
          attr(value, "row.names") <- rn2
        }

    value
  }

###### Other optmatch.dlist common methods #####
#' @export
dim.optmatch.dlist <- function(x) {
  dims <- lapply(x, dim)
  return(Reduce(function(x,y) { c(x[1] + y[1], x[2] + y[2])}, dims, c(0,0)))
}

#' @export
dimnames.optmatch.dlist <- function(x) {
  dnms <- lapply(x, dimnames)
  return(Reduce(function(x,y) {list(treated = c(x$treated, y[[1]]), control =
  c(x$control, y[[2]]))}, dnms, list(treated = c(), control = c())))
}

#' @export
as.matrix.optmatch.dlist <- function(x, ...) {
  xdim <- dim(x)
  tmp <- matrix(Inf, nrow = xdim[1], ncol = xdim[2], dimnames = dimnames(x))

  for (i in 1:length(x)) {
    submatrix <- x[[i]]
    subrows <- rownames(submatrix)
    subcols <- colnames(submatrix)
    tmp[subrows, subcols] <- submatrix
  }
  return(tmp)
}

#' @export
subset.optmatch.dlist <- function(x, subset, select, ...) {
  subset(as.matrix(x), subset, select, ...)
}


#' @rdname subdim-methods
#' @export
subdim.optmatch.dlist <- function(x) {
  lapply(x, dim)
}
