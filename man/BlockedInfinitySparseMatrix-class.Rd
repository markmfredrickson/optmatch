\name{BlockedInfinitySparseMatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{BlockedInfinitySparseMatrix-class}

\title{Class \code{"BlockedInfinitySparseMatrix"}}
\description{
  Like the parent class \code{"InfinitySparseMatrix"}, the
  \code{"BlockedInfinitySparseMatrix"} encodes a sparse matching problem. This
  class extends the parent with completely disconnected strata. This allows for
  splitting up large matching problems into smaller subproblems that can each
  be solved more quickly and efficiently.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("BlockedInfinitySparseMatrix", ...)}, but users would probably prefer
  to use \code{\link{exactMatch}}, a convenience function to create stratified
  sparse matching problems.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}: The real valued entries of the matrix. Default data when accesing the object.}
    \item{\code{groups}:}{Object of class \code{"factor"}: A named factor indicating the group memebership of the treated and control units. }
    \item{\code{cols}:}{Object of class \code{"numeric"}: The column indices of real entries. }
    \item{\code{rows}:}{Object of class \code{"numeric"}: The row indices of the real entries. }
    \item{\code{dimension}:}{Object of class \code{"numeric"}: A two element vector of height and width. }
    \item{\code{colnames}:}{Object of class \code{"OptionalCharacter"} Names of the columns. }
    \item{\code{rownames}:}{Object of class \code{"OptionalCharacter"} Names of the rows. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{InfinitySparseMatrix}"}, directly.
Class \code{"\linkS4class{numeric}"}, by class "InfinitySparseMatrix", distance 2.
Class \code{"\linkS4class{DistanceSpecification}"}, by class "InfinitySparseMatrix", distance 2.
Class \code{"\linkS4class{vector}"}, by class "InfinitySparseMatrix", distance 3.
}
\section{Methods}{
No public methods are defined specifically for class "BlockedInfinitySparseMatrix". All of the methods in \code{linkS4class{InfinitySparseMatrix}} work for this class as well.
}
\author{Mark M. Fredrickson}
\note{The levels in the \code{groups} factor will later show up as the prefixes in the results of \code{\link{fullmatch}}.}

\seealso{
\code{\link{exactMatch}}, \code{\linkS4class{InfinitySparseMatrix}} for links to other classes ~~~
}
\keyword{classes}
