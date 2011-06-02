\name{makedist}
\alias{makedist}
\title{Assemble match distances from a data frame}
\description{
  Helper function to produce first arguments to \code{fullmatch()}, reducing memory requirements
  for \code{fullmatch()} and heading off certain user errors.
}
\usage{
makedist(structure.fmla, data, fn = function(trtvar, dat, ...) {
    matrix(0, sum(trtvar), sum(!trtvar), dimnames = list(names(trtvar)[trtvar], names(trtvar)[!trtvar]))
}, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure.fmla}{A formula defined w.r.t. data frame \code{data},
  with a treatment variable on the LHS and either \sQuote{1}, if no
  stratification prior to matching, or variables
defining pre-matching stratification on the RHS.}
  \item{data}{ A data frame in which \code{structure.fmla} is evaluated
    and \code{fn} operates. }
  \item{fn}{ A user-supplied function to compute distances. See
    details and examples.}
  \item{\dots}{ Additional arguments to \code{fn} }
}
\details{
  \code{fn} should be a function with first two arguments
  \code{trtvar}, a treatment variable, and 
  \code{dat}, a data frame.  There may be additional arguments.  If the
  function uses variables in \code{dat}, these should be referenced
  using names from the input \code{trtvar}, particularly if the sample is being
  split into strata (ie \code{structure.fmla} has a non-trivial RHS).
  When this happens, \code{fn} will be passed a \code{trtvar} input
  observations for only a subset of the rows of \code{dat}, so it has to
  use \code{trtvar} to decide which rows of \code{dat} to operate on; it
  does this by lining up names of the (shorter) vector \code{trtvar}
  with row names of \code{dat}.
}
\value{
  A list of matrices, one for each subclass defined by the interaction
  of variables appearing on the right hand side of
  \code{structure.fmla}.  Each of these is a number of treatments by
  number of controls matrix of distances, with the distance between
  treatments and controls calculated by the user-given function \code{fn}.

  The list also has some attributes that are not of direct interest to
  the user, but are used by \code{fullmatch()}.
}
\author{Ben Hansen }
\note{ Use of this function to prepare input to \code{fullmatch()}
  prevents two common problems.  First, the function encourages you to
  stratify a large data set and match within strata, then calculates
  distances only between potential matches within the same stratum.
  Whether or not you stratify, the function has the advantage of keeping
track of the order of observations in a data frame from which distances
are generated.  It passes this info to \code{fullmatch()}, which makes
sure to return the vector a vector with the ordering of the generating
data frame.  It's still possible to create matrices of distances without
this function, but then \code{fullmatch()} has no way of knowing the
order of observations in whatever data frame you're working from, since
that info is lost in the transition to a distance matrix.}

\seealso{\code{\link{fullmatch}}, \code{\link{pscore.dist}}, \code{\link{mahal.dist}}}
\examples{
data(nuclearplants)

##-- A distance function used in P. Rosenbaum's (2002) book
rankdiffs <- function(trtvar, dat, vars)
{
dmt <- matrix(0,sum(trtvar), sum(!trtvar)) 
for (vv in vars) {
vvr <- rank(dat[names(trtvar),vv])
dmt <- dmt + abs(outer(vvr[trtvar], vvr[!trtvar],"-"))
}
round(dmt)              
}
##-- Gives a warning because this fn doesn't assign dimnames
(rdd1 <- makedist(pr~1, nuclearplants[nuclearplants$pt==0,], rankdiffs, c("cap","date")))
fullmatch(rdd1)
##-- fullmatch() knows its value should be ordered as the nuclearplants data set is
rdd1$m
##-- now fullmatch() doesn't know the proper order of units and has to guess
fullmatch(rdd1$m)
(rdd2 <- makedist(pr~pt, nuclearplants, rankdiffs, c("cap","date")))
fullmatch(rdd2)

##- Distance on a propensity score
scalardiffs <- function(trtvar,data,scalarname) {
sclr <- data[names(trtvar), scalarname]
names(sclr) <- names(trtvar)
abs(outer(sclr[trtvar],sclr[!trtvar], '-'))
}
nuclearplants$pscore <- glm(pr~.-(pr+cost), family=binomial(),
                      data=nuclearplants)$linear.predictors
##-- Distance for propensity score matching w/o prior stratification
psd1 <- makedist(pr~1, nuclearplants, scalardiffs, "pscore")
fullmatch(psd1)
##-- Distance for propensity score matching within levels of "pt"
psd2 <- makedist(pr~pt, nuclearplants, scalardiffs, "pscore")
fullmatch(psd2)

}
\keyword{nonparametric }% at least one, from doc/KEYWORDS

