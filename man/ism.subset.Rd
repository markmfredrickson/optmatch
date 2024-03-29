% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InfinitySparseMatrix.R
\name{subset.InfinitySparseMatrix}
\alias{subset.InfinitySparseMatrix}
\alias{[,InfinitySparseMatrix-method}
\alias{[<-,InfinitySparseMatrix-method}
\title{Subsetting for InfinitySparseMatrices}
\usage{
\method{subset}{InfinitySparseMatrix}(x, subset, select, ...)

\S4method{[}{InfinitySparseMatrix}(x, i, j = NULL, ..., drop = TRUE)

\S4method{[}{InfinitySparseMatrix}(x, i, j) <- value
}
\arguments{
\item{x}{InfinitySparseMatrix to be subset or bound.}

\item{subset}{Logical expression indicating rows to keep.}

\item{select}{Logical expression indicating columns to keep.}

\item{...}{Other arguments are ignored.}

\item{i}{Row indices.}

\item{j}{Col indices.}

\item{drop}{Ignored.}

\item{value}{replacement values}
}
\value{
An InfinitySparseMatrix with only the selected elements.
}
\description{
This matches the syntax and semantics of
subset for matrices.
}
\author{
Mark Fredrickson
}
