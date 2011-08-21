\name{caliper}
\alias{caliper}
\title{Prepare matching distances suitable for matching within calipers}
\description{
  Encodes calipers, or maximum allowable distances within which to
  match, calling \code{mdist()} on the supplied arguments to calculate
  the distance on which the caliper is based. (One then applies
  \code{pairmatch()} or \code{fullmatch()} to the value of
  \code{caliper}, or to the sum of it and another distance, to match
  within calipers.)
}
\usage{
caliper(width, ..., exclude = c(), compare = `<=`)
}
\arguments{
  \item{width}{The width of the caliper: how wide of a margin to allow in matches.}
  \item{...}{Arguments are passed to \code{mdist()}, which forms the basic distance matrix. The caliper is then fit to this matrix. See \code{\link{mdist}} for more information.}
  \item{exclude}{A character vector of observations (corresponding to rownames) to exclude from the caliper.}
  \item{cmpare}{A function that decides that whether two observations are with the caliper. The default is \code{`<=`};
  you may prefer \code{`<`}}
}
\details{
    
  \code{width} provides the size of the caliper, the allowable distance for
  matching. The unit of width will depend on the additional arguments, used 
  to create a basic distance matrix using \code{\link{mdist}}. For a
  forumla, width will be in the unit of the selected covariate. For a GLM
  \code{x} argument, the width will be relative to the computed propensity
  score. Careful consideration should be given to what this unit means in
  your analysis. Similarly, a previous matching matrix will use those units.

}
\value{
  Object of class \code{DistanceSpecification}, which is suitable to be given
  as \code{distance} argument to \code{\link{fullmatch}} or
  \code{\link{pairmatch}}. The caliper will be only zeros and \code{Inf} values,
  indicating a possible match or no possible match, respectively.

}
\author{Mark M. Fredrickson and Ben B. Hansen}

\references{
  P.~R. Rosenbaum and D.~B. Rubin (1985),
\sQuote{Constructing a control group using multivariate matched sampling
  methods that incorporate the propensity score},
 \emph{The American Statistician}, \bold{39} 33--38.
}

\seealso{\code{\link{exactMatch}}, \code{\link{mdist}}, \code{\link{fullmatch}}, \code{\link{pairmatch}}}

\examples{
data(nuclearplants)

### Caliper  of .2 pooled SDs in the propensity score
ppty <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
(pptycaliper <- caliper(ppty, width = .2))

identical(pptycaliper,  # If not writing 'width=',
caliper(.2,ppty))        # give your width first.

### Caliper on a pre-formed distance
ppty.dist <- mdist(ppty) 
identical(caliper(ppty.dist, width = .2), pptycaliper)

### caliper on the Mahalanobis distance
caliper(width = 3, pr ~ t1 + t2, data = nuclearplants)             

### Mahalanobis distance matching with a caliper
### of 1 pooled SD in the propensity score:
( mhd.pptyc <- caliper(ppty, width = 1) +
  mdist(pr ~ t1 + t2, data = nuclearplants) )
pairmatch(mhd.pptyc)

### Excluding observations from caliper requirements:
caliper(width = 3, pr ~ t1 + t2, data = nuclearplants, exclude = c("A", "f"))

}
\keyword{nonparametric }% at least one, from doc/KEYWORDS
