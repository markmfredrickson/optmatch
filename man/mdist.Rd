\name{mdist}
\alias{mdist}
\alias{mdist-methods}
\alias{mdist,bigglm-method}
\alias{mdist,formula-method}
\alias{mdist,function-method}
\alias{mdist,glm-method}
\alias{mdist,InfinitySparseMatrix-method}
\alias{mdist,matrix-method}
\alias{mdist,numeric-method}
\title{Create matching distances}
\description{
  A generic function, with several supplied methods, for creating
  matrices of distances between observations to be used in the match process.
  Using these matrices, 
  \code{pairmatch()} or \code{fullmatch()} can determine optimal matches.
}
\usage{
mdist(x, exclusions = NULL, ...)

\S4method{mdist}{function}(x, exclusions = NULL, z = NULL, data = NULL, ...) 

\S4method{mdist}{formula}(x, exclusions = NULL, data = NULL, subset = NULL, inv.scale.matrix = NULL, COV = cov, ...)

\S4method{mdist}{glm}(x, exclusions = NULL, standardization.scale = mad, ...)

\S4method{mdist}{bigglm}(x, exclusions = NULL, data = NULL, standardization.scale = mad, ...)

\S4method{mdist}{numeric}(x, exclusions = NULL, ...)

\S4method{mdist}{InfinitySparseMatrix}(x, exclusions = NULL, ...)

\S4method{mdist}{matrix}(x, exclusions = NULL, ...)
}
\arguments{
  \item{x}{The object to use as the basis for dispatch to the specialized methods}
  \item{exclusions}{An optional distance specification (perhaps created by \code{exactMatch} or \code{caliper}) that defines for which entries values should be computed. This can significantly speed up computation for sparse matching problems.}
  \item{z}{A vector indicating the whether the entries in \code{data} are in the treated or control groups.}
  \item{data}{A \code{data.frame} or \code{matrix} containing variables used by the method to construct the distance matrix.}
  \item{subset}{A subset of the data to use in creating the distance specification.}
  \item{inv.scale.matrix}{A matrix that scales the distance computation. The default Mahalanobis distance scales the squared Euclidean distance by
the inverse of the covariance matrix. Other scale matrices can be supplied, for example,
the identity matrix will result in squared Euclidean distance.}
  \item{COV}{A covariance computing function that can be supplied to the \code{formula} method.}
  \item{standardization.scale}{Standardizes the data based on the median absolute deviation (by default).}
  \item{...}{Additional method arguments.}
}
\details{
  The \code{function} method takes as its \code{x} argument a function of two arguments.
  This function should in turn expect two arguments, a \code{data.frame} carrying information
  about treatment group members and a \code{data.frame} carrying information about control units. 
  The \code{function} method also expects a \code{z} argument, a vector indicating whether each unit is in the treated
  or control groups; on the basis of this vector the \code{data} argument is what is chopped up and fed
  to the function passed to \code{mdist} as argument \code{x}. While this might sound
  complicated at first, it's quite flexible and, once you're used to it, easy to use. For example, the simple
  1-dimensional distance of treated and control units can be implemented as:
  \code{d <- 1:20 ; names(d) <- letters[1:20] ; mdist(`-`, z = rep(c(T,F),10), data = d)}.  See also additional examples below.  

  The formula method produces, by default, a Mahalanobis distance specification
  based on the formula \code{treatment ~ var1 + var2 + ... }, where
  \code{treatment} iexclusions = NULLs the treatment indicator. A Mahalanobis
  distance scales the squared Euclidean distance by the inverse of the
  covariance matrix. Other scale matrices can be supplied, for example, the
  identity matrix will result in squared Euclidean distance. If you wish to use
  an alternative function to compute the covariance, pass it in the \code{COV}
  argument.

  An \code{glm} method takes an argument of class \code{glm} as
  first argument.  It assumes that this object is a fitted propensity
  model, extracting distances on the linear propensity score (logits of
  the estimated conditional probabilities) and, by default, rescaling the distances
  by the reciprocal of the pooled s.d. of treatment- and control-group propensity scores.
  (The scaling uses \code{mad}, for resistance to outliers, by default; this can be
  changed to the actual s.d., or rescaling can be skipped entirely, by
  setting argument \code{standardization.scale} to \code{sd} or
  \code{NULL}, respectively.)  The \code{bigglm}
  method works analogously with \code{bigglm} objects, created by
  the \code{bigglm} function from package \sQuote{biglm}, which can
  handle bigger data sets than the ordinary glm function can. 

  The \code{numeric} method is simply a placeholder that returns an error
  message and some helpful suggestions on how to create distance matrices, as
  there is no obvious way to create distance matrix from just a numeric vector.

  The \code{matrix} and \code{InfinitySparseMatrix} just return their arguments as these objects are already valid distance specifications.
}
\value{
  Object of the class union \code{DistanceSpecification} (either a
  \code{matrix} or a \code{InfinitySparseMatrix} or derived class), which is
  suitable to be given as \code{distance} argument to \code{\link{fullmatch}}
  or \code{\link{pairmatch}}. 

}
\author{ Mark M. Fredrickson} %Ben B. Hansen

\references{
  P.~R. Rosenbaum and D.~B. Rubin (1985),
\sQuote{Constructing a control group using multivariate matched sampling
  methods that incorporate the propensity score},
 \emph{The American Statistician}, \bold{39} 33--38.
}

\seealso{\code{\link{fullmatch}}, \code{\link{pairmatch}}, \code{\link{caliper}}, \code{\link{exactMatch}}}
\examples{
data(nuclearplants)
mdist.examples <- list()
### Propensity score distances.
### Recommended approach:
(aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants))
mdist.examples$ps1 <- mdist(aGlm)
### A second approach: first extract propensity scores, then separately
### create a distance from them.  (Useful when importing propensity
### scores from an external program.)
plantsPS <- predict(aGlm)
mdist.examples$ps2 <- mdist(pr~plantsPS, data=nuclearplants)^(1/2)
### Full matching on the propensity score.
fullmatch(mdist.examples$ps1)
fullmatch(mdist.examples$ps2)
### Because mdist.glm uses robust estimates of spread, 
### the results differ in detail -- but they are close enough
### to yield similar optimal matches.
all(fullmatch(mdist.examples$ps1)==fullmatch(mdist.examples$ps2)) # The same

### Mahalanobis distance:
mdist.examples$mh1 <- mdist(pr ~ t1 + t2, data = nuclearplants)

### Absolute differences on a scalar:

sdiffs <- function(treatments, controls) {
  abs(treatments$t1 - controls$t1)
}

(absdist <- mdist(sdiffs, z = nuclearplants$pr, data = nuclearplants, 
                  exclusions = exactMatch(pr ~ pt, nuclearplants)))

### Pair matching on the variable `t1`:
pairmatch(absdist)


### Propensity score matching within subgroups:
mdist.examples$ps3 <- mdist(aGlm, exactMatch(pr ~ pt, nuclearplants))
fullmatch(mdist.examples$ps3)

### Propensity score matching with a propensity score caliper:
mdist.examples$pscal <- mdist.examples$ps1 + caliper(mdist.examples$ps1, 1)
fullmatch(mdist.examples$pscal) # Note that the caliper excludes some units

### A Mahalanobis distance for matching within subgroups:
mdist.examples$mh2 <- mdist(pr ~ t1 + t2 , data = nuclearplants,
                            excludes = exactMatch(pr ~ pt, nuclearplants))

### Mahalanobis matching within subgroups, with a propensity score
### caliper:
fullmatch(mdist.examples$mh2 + caliper(mdist.examples$ps3, 1))

}

\keyword{nonparametric }% at least one, from doc/KEYWORDS
