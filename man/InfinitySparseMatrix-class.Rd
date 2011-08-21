\name{InfinitySparseMatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{InfinitySparseMatrix-class}
\alias{Arith,InfinitySparseMatrix,BlockedInfinitySparseMatrix-method}
\alias{Arith,InfinitySparseMatrix,InfinitySparseMatrix-method}
\alias{Arith,InfinitySparseMatrix,matrix-method}
\alias{Arith,matrix,InfinitySparseMatrix-method}
\alias{as.matrix,InfinitySparseMatrix-method}
\alias{coerce,InfinitySparseMatrix,matrix-method}
\alias{coerce,matrix,InfinitySparseMatrix-method}
\alias{colnames,InfinitySparseMatrix-method}
\alias{colnames<-,InfinitySparseMatrix-method}
\alias{dimnames,InfinitySparseMatrix-method}
\alias{dimnames<-,InfinitySparseMatrix-method}
\alias{mdist,InfinitySparseMatrix-method}
\alias{prepareMatching,InfinitySparseMatrix-method}
\alias{rownames,InfinitySparseMatrix-method}
\alias{rownames<-,InfinitySparseMatrix-method}

\title{Class \code{"InfinitySparseMatrix"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("InfinitySparseMatrix", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cols}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rows}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dimension}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{colnames}:}{Object of class \code{"OptionalCharacter"} ~~ }
    \item{\code{rownames}:}{Object of class \code{"OptionalCharacter"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{numeric}"}, from data part.
Class \code{"\linkS4class{DistanceSpecification}"}, directly.
Class \code{"\linkS4class{vector}"}, by class "numeric", distance 2.
}
\section{Methods}{
  \describe{
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "BlockedInfinitySparseMatrix")}: ... }
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix")}: ... }
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "matrix")}: ... }
    \item{Arith}{\code{signature(e1 = "matrix", e2 = "InfinitySparseMatrix")}: ... }
    \item{as.matrix}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "InfinitySparseMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "InfinitySparseMatrix")}: ... }
    \item{colnames}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{colnames<-}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{dimnames}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{dimnames<-}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{mdist}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{prepareMatching}{\code{signature(distances = "InfinitySparseMatrix")}: ... }
    \item{rownames}{\code{signature(x = "InfinitySparseMatrix")}: ... }
    \item{rownames<-}{\code{signature(x = "InfinitySparseMatrix")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("InfinitySparseMatrix")
}
\keyword{classes}
