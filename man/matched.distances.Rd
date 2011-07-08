\name{matched.distances}
\alias{matched.distances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine distances between matched units }
\description{
  From a match (as produced by \code{pairmatch} or \code{fullmatch}) and
  a distance, extract the distances of matched units from their matched counterparts.
}
\usage{
matched.distances(matchobj, distance,preserve.unit.names=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matchobj}{ Value of a call to \code{pairmatch} or \code{fullmatch}. }
  \item{distance}{ Either a distance matrix or the value of a call to
    or \code{makedist}.}
  \item{preserve.unit.names}{Logical.  If true, for each matched set
    \code{matched.distances} returns the submatrix of the distance
    matrix corresponding to it; if false, a vector containing the
    distances in that submatrix is returned.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of numeric vectors (or matrices) of distances, one for each matched set.
  Note that a matched set with 1 treatment and k controls, or with k
  treatments and 1 control, has k, not k+1, distances.
}
\author{ Ben B. Hansen }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\examples{
data(plantdist)
plantsfm <- fullmatch(plantdist) 
(plantsfm.d <- matched.distances(plantsfm,plantdist,pres=TRUE))
unlist(lapply(plantsfm.d, max))
mean(unlist(plantsfm.d))
}
\keyword{nonparametric }% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
