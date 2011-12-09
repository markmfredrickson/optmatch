\name{minControlsCap}
\alias{minControlsCap}
\alias{maxControlsCap}

\title{Set thinning and thickening caps for full matching}
\description{
  Functions to find the largest value of min.controls,
  or the smallest value of max.controls, for which a full matching
  problem is feasible.  These are determined by
  constraints embedded in the matching problem's distance matrix.
}
\usage{
minControlsCap(distance, max.controls = NULL)
maxControlsCap(distance, min.controls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{ Either a matrix of nonegative, numeric discrepancies,
  or a list of such matrices. (See \code{\link{fullmatch}} for details.)}
  \item{max.controls}{ Optionally, set limits on the maximum number of
    controls per matched set.  (Only makes sense for \code{minControlsCap}.)}
  \item{min.controls}{ Optionally, set limits on the minimum number of
    controls per matched set.  (Only makes sense for \code{maxControlsCap}.)}
}
\details{
    The function works by repeated application of full matching, so on
  large problems it can be time-consuming.
}
\value{
  For \code{minControlsCap}, \code{strictest.feasible.min.controls} and
\code{given.max.controls}. For \code{maxControlsCap},
\code{given.min.controls} and \code{strictest.feasible.max.controls}.

  \item{strictest.feasible.min.controls}{The largest values of
    the \code{\link{fullmatch}} argument \code{min.controls} that yield
    a full match;}
  \item{given.max.controls}{
    The \code{max.controls} argument given to \code{minControlsCap}
    or, if none was given, a vector of \code{Inf}s. }
  \item{given.min.controls}{
    The \code{min.controls} argument given to \code{maxControlsCap}
    or, if none was given, a vector of \code{0}s; }
  \item{strictest.feasible.max.controls}{The smallest values of
    the \code{\link{fullmatch}} argument \code{max.controls} that yield
    a full match.}
}
\references{Hansen, B.B. and S. Olsen Klopfer (2006),
  \sQuote{Optimal full matching and related designs via network flows},
  \emph{Journal of Computational and Graphical Statistics} \bold{15}, 609--627.  }
\author{Ben B. Hansen }
\note{ Essentially this is just a line search.  I've done several things
  to speed it up, but not everything that might be done.  At present,
  not very thoroughly tested either: you might check the final results
    to make sure that \code{\link{fullmatch}} works with the values of
    \code{min.controls} (or \code{max.controls}) suggested by these
    functions, and that it ceases to work if you increase (decrease)
    those values.  Comments appreciated.}


\seealso{\code{\link{fullmatch}}}
\examples{
### TEMPORARY COMMENT OUT BY MMF
### ##---- Should be DIRECTLY executable !! ----
### ##-- ==>  Define data, use random,
### ##--	or do  help(data=index)  for the standard data sets.
### plantdist <- matrix(nrow=7, ncol=19,byrow=TRUE,data=c(
### 28, 0, 3,22,14,30,17,28,26,28,20,22,23,26,21,18,34,40,28,
### 24, 3, 0,22,10,27,14,26,24,24,16,19,20,23,18,16,31,37,25,
### 10,18,14,18, 4,12, 6,11, 9,10,14,12, 6,14,22,10,16,22,28,
###  7,28,24, 8,14, 2,10, 6,12, 0,24,22, 4,24,32,20,18,16,38,
### 17,20,16,32,18,26,20,18,12,24, 0, 2,20, 6, 8, 4,14,20,14,
### 20,31,28,35,20,29,22,20,14,26,12, 9,22, 5,15,12, 9,11,12,
### 14,32,29,30,18,24,17,16,10,22,12,10,17, 6,16,14, 4, 8,17),
### dimnames=list(c("A","B","C","D","E","F","G"),
### c("H","I","J","K","L","M","N","O","P","Q","R",
### "S","T","U","V","W","X","Y","Z")))
### 
### (tmn <- minControlsCap(plantdist)$strictest)
### maxControlsCap(plantdist, min=tmn)
### ### START HERE! Need to recreate split dist from the plantdist matrix above.
### splitdist <- matrix(Inf, nrow = 7, ncol = 19)
### splitdist[1:3, 1:9] <- 0
### splitdist[4:7, 10:19] <- 0
### (tmn <- minControlsCap(splitdist + plantdist)$strictest)
### maxControlsCap(splitdist + plantdist, min=tmn)
### }
\keyword{ optimize }% at least one, from doc/KEYWORDS
