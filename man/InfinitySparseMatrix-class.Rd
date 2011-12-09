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
\alias{prepareMatching,InfinitySparseMatrix-method}
\alias{rownames,InfinitySparseMatrix-method}
\alias{rownames<-,InfinitySparseMatrix-method}
\alias{cbind.InfinitySparseMatrix}
\alias{rbind.InfinitySparseMatrix}
\alias{dim.InfinitySparseMatrix}
\alias{subset.InfinitySparseMatrix}
\alias{t.InfinitySparseMatrix}

\title{Class \code{"InfinitySparseMatrix"}}
\description{
The \code{InfinitySparseMatrix} class is a specialized matrix class where the default
(and most probably most common) entry is Infinity rather than zero. This matrix is available
as a \code{DistributionSpecification} object to specify sparse matching problems. Coercion to 
dense matrices, element wise arithmetic, column and row binding methods exist.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("InfinitySparseMatrix", ...)};
however, users will prefer the helper function \code{\link{makeInfinitySparseMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}: The real valued entries of the matrix. Default data when accesing the object.}
    \item{\code{cols}:}{Object of class \code{"numeric"}: The column indices of real entries. }
    \item{\code{rows}:}{Object of class \code{"numeric"}: The row indices of the real entries. }
    \item{\code{dimension}:}{Object of class \code{"numeric"}: A two element vector of height and width. }
    \item{\code{colnames}:}{Object of class \code{"OptionalCharacter"}: Names of the columns. }
    \item{\code{rownames}:}{Object of class \code{"OptionalCharacter"}: Names of the rows. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{numeric}"}, from data part. Most numeric functions are valid.
Class \code{"\linkS4class{DistanceSpecification}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "BlockedInfinitySparseMatrix")}: Element wise arithmetic. }
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix")}: Element wise arithmetic. }
    \item{Arith}{\code{signature(e1 = "InfinitySparseMatrix", e2 = "matrix")}: Element wise artihmetic, infinite entries removed. }
    \item{Arith}{\code{signature(e1 = "matrix", e2 = "InfinitySparseMatrix")}: Element wise artihmetic, infinite entries removed. }
    \item{as.matrix}{\code{signature(x = "InfinitySparseMatrix")}: Convert to a dense matrix. }
    \item{coerce}{\code{signature(from = "InfinitySparseMatrix", to = "matrix")}: Convert to a dense matrix. }
    \item{coerce}{\code{signature(from = "matrix", to = "InfinitySparseMatrix")}: Convert a dense matrix to an InfinitySparseMatrix. }
    \item{colnames}{\code{signature(x = "InfinitySparseMatrix")}: Get column names. }
    \item{colnames<-}{\code{signature(x = "InfinitySparseMatrix")}: Set column names. }
    \item{dimnames}{\code{signature(x = "InfinitySparseMatrix")}: Get row and column names. }
    \item{dimnames<-}{\code{signature(x = "InfinitySparseMatrix")}: Set row and column names. }
    \item{prepareMatching}{\code{signature(distances = "InfinitySparseMatrix")}: Turn a InfinitySparseMatrix into a cannonical matching problem. }
    \item{rownames}{\code{signature(x = "InfinitySparseMatrix")}: Get the row names. }
    \item{rownames<-}{\code{signature(x = "InfinitySparseMatrix")}: Set the row names. }
    \item{cbind.InfinitySparseMatrix}{\code{signature(x = "InfinitySparseMatrix", y = "DistanceSpecification")}: Column bind. }
    \item{rbind.InfinitySparseMatrix}{\code{signature(x = "InfinitySparseMatrix", y = "DistanceSpecification")}: Row bind. }
    \item{dim.InfinitySparseMatrix}{\code{signature(x = "InfinitySparseMatrix")}: Get row and column names. }
    \item{subset.InfinitySparseMatrix}{\code{signature(x = "InfinitySparseMatrix", subset, ...)}: Create an InfinitySparseMatrix of the same dimenions but only containing the entries specified in \code{subset}. }
    \item{t.InfinitySparseMatrix}{\code{signature(x = "InfinitySparseMatrix")}: Matrix transpose. }

	 }
}
% \references{
%%  ~~put references to the literature/web site here~~
% }
\author{Mark M. Fredrickson}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
\code{\linkS4class{DistanceSpecification}}
}
\examples{
showClass("InfinitySparseMatrix")
}
\keyword{classes}
