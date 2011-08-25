\name{makedist}
\alias{makedist}
\title{Assemble match distances from a data frame}
\description{
  Helper function to produce first arguments to \code{fullmatch()}, reducing memory requirements
  for \code{fullmatch()} and heading off certain user errors. You probably want to use \code{mdist}
  instead.
}
\usage{
makedist(z, data, distancefn, exclusions = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{A treatment indicator variable with precisely two unique values}
  \item{data}{A vector, data frame or matrix the same number of observations as \code{z}}
  \item{distancefn}{ A user-supplied function to compute distances.
    takes two arguments: the treated data for each treatment-control pair and the control data for each pair.
    Returns a list of distances the same length of the pairs objects.
  \item{exclusions}{ A (sparse) distance specification that indicates which pairs are valid matches. }
}
\details{
}
\value{
  A DistanceSpecification with the computed distances between treated and control units.
}
\author{Ben Hansen and Mark M. Fredrickson}
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

\seealso{\code{\link{fullmatch}}, \code{\link{mdist}} }
\examples{
data(nuclearplants)

##-- A distance function used in P. Rosenbaum's (2002) book
ranked.plants <- apply(nuclearplants, 2, rank)

rankdiffs <- function(treateds, controls) {
  # treateds and controls are data.frames of equal height
  rowSums(abs(treateds - controls))
}

(rdd <- makedist(nuclearplants$pr, ranked.plants[, c("cap", "date")], rankdiffs))
fullmatch(rdd)

##- Distance on a propensity score
scalar.diffs <- function(treated, controls) {
  # here treated and controls are two equal length vectors
  abs(treated - controls)
}

pscore <- glm(pr~.-(pr+cost), family=binomial(),
                      data=nuclearplants)$linear.predictors

names(pscore) <- rownames(nuclearplants)

##-- Distance for propensity score matching w/o prior stratification
psd1 <- makedist(nuclearplants$pr, pscore, scalar.diffs)
fullmatch(psd1)

##-- Distance for propensity score matching within levels of "pt"
psd2 <- makedist(nuclearplants$pr, pscore, scalar.diffs,
                 exactMatch(pr ~ pt, data = nuclearplants))
fullmatch(psd2)

}
\keyword{nonparametric }% at least one, from doc/KEYWORDS

