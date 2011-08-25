\name{exactMatch-methods}
\docType{methods}
\alias{exactMatch-methods}
\alias{exactMatch,formula-method}
\alias{exactMatch,vector-method}
\title{Methods for Function \code{exactMatch} in Package \pkg{optmatch}}
\description{
  An exact match is one based on a factor. Within a level, all observations are allowed to be matched.
  An exact match can be combined with another distance matrix to create a set of matching subproblems.

}

\section{Methods}{
\describe{

\item{\code{signature(x = "formula", data = NULL, subset = NULL, na.action = NULL, ...)}}{
  The formula should be of the form \code{treatment ~ factor1 + factor2 + ...}. The optional arguments are passed to \code{model.frame} before creating the interation of all factors.
} 

\item{\code{signature(x = "vector", treatment)}}{
  In addition to a vector \code{x} indicating the factor level, you must all include a second vector of equal length \code{treatment} which indicates which units are treatment and which are control.
}
}}
\keyword{methods}
