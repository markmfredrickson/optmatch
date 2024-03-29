% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InfinitySparseMatrix.R
\name{sort.InfinitySparseMatrix}
\alias{sort.InfinitySparseMatrix}
\alias{sort.BlockedInfinitySparseMatrix}
\title{Sort the internal structure of an InfinitySparseMatrix.}
\usage{
\method{sort}{InfinitySparseMatrix}(x, decreasing = FALSE, ..., byCol = FALSE)

\method{sort}{BlockedInfinitySparseMatrix}(x, decreasing = FALSE, ..., byCol = FALSE)
}
\arguments{
\item{x}{An \code{InfinitySparseMatrix} or
\code{BlockedInfinitySparseMatrix}.}

\item{decreasing}{Logical. Should the sort be increasing or decreasing?
Default \code{FALSE}.}

\item{...}{Additional arguments ignored.}

\item{byCol}{Logical. Defaults to \code{FALSE}, so the returned ISM is
row-dominant. \code{TRUE} returns a column-dominant ISM.}
}
\value{
An object of the same class as \code{x} which is sorted according to
  \code{byCol}.
}
\description{
Internally, an \code{InfinitySparseMatrix} (Blocked or non) comprises of
vectors of values, row positions, and column positions. The ordering of
these vectors is not enforced. This function sorts the internal structure,
leaving the external structure unchanged (e.g. \code{as.matrix(ism)} and
\code{as.matrix(sort(ism))} will look identical despite sorting.)
}
\details{
By default, the \code{InfinitySparseMatrix} is row-dominant, meaning the row
positions are sorted first, then column positions are sorted within each
row. Use argument \code{byCol} to change this.
}
