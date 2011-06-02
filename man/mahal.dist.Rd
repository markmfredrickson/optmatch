\name{mahal.dist} 
\alias{mahal.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Assemble Mahalanobis distances and prepare for matching }
\description{
  Calculates squared Mahalanobis distances between treatment and control
  observations on given variables,
  assembling them into a discrepancy matrix (or matrices) from which
  \code{pairmatch()} or \code{fullmatch()} can determine optimal
  matches.  (If vectors \code{x} and \code{y} encode two observations'
  values of the specified variables, then the squared Mahalanobis
  distance between them is
  \deqn{
    D^2 = (x - y)'\Sigma^{-1}(x-y)
    }
  .)
}
\usage{
mahal.dist(distance.fmla, data, structure.fmla = NULL, inverse.cov = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance.fmla}{A formula with variables to be combined in
    the Mahalanobis distance on its right-hand side and the
  treatment variable on its left.}
  \item{data}{ Data frame in which \code{distance.fmla} and (if
    specified) \code{structure.fmla} are to be evaluated. }
  \item{structure.fmla}{ Optional formula argument specifying subclasses
  within which matches are to be performed.  If omitted, no
  subclassification is done. If it is given, its left-hand side gives
  the treatment variable and its RHS  gives
    variables on which to stratify the sample prior to matching. }
  \item{inverse.cov}{ The inverse covariance of the variables to be
    combined into the Mahalanobis distance (optional). }
}
\details{
   Mahalanobis distance tracks the discrepancy between points on a
   number of given variables, after standardizing the variables and
   taking account of their covariances.  It is best suited to variables
   whose joint distribution resembles a multivariate Normal.

     The purpose of giving a \code{structure.fmla} argument is to speed
  up large problems.  Variables appearing on its
  right-hand side will be interacted to create the subclasses.  If
  \code{structure.fmla} is given then its LHS is used to define
  treatment and control groups (and one doesn't have to put anything on
  the LHS of \code{distance.fmla}).

  The function attempts to calculate the inverse covariance itself, so
  ordinarily you shouldn't need to give it one.  If you'll be calling
  the function repeatedly, however, it may speed things up to compute
  and store the inverse covariance once, rather than each time this
  function is called; in that case you can save time by giving the
  \code{inverse.covariance} argument.
}
\value{
  Object of class \code{optmatch.dlist}, which is suitable to be given
  as \code{distance} argument to \code{\link{fullmatch}} or
  \code{\link{pairmatch}}.

  Specifically, a list of matrices, one for each subclass defined by the
  interaction of variables appearing on the right hand side of
  \code{structure.fmla}.  Each of these is a number of treatments by
  number of controls matrix of propensity distances.
  The list also carries some metadata as attributes, data that is not of direct interest to
  the user but is useful to \code{fullmatch()} and \code{pairmatch()}.
}

\author{ Ben B. Hansen}

\seealso{ \code{\link{makedist}}, \code{\link{pscore.dist}}, \code{\link{fullmatch}}, \code{\link{pairmatch}} }

\examples{
data(nuclearplants)
mhd1 <- mahal.dist(pr~date+cum.n, nuclearplants)
lapply(mhd1, round)
attributes(mhd1)
fullmatch(mhd1)
##- Mahalanobis within subclasses defined by levels of pt
mhd2 <- mahal.dist(~date+cum.n, nuclearplants, pr~pt)
lapply(mhd2, round)
fullmatch(mhd2)
##- Trick mahal.dist into returning absolute differences on a scalar.
mhd3 <- mahal.dist(pr~date, nuclearplants,
inverse.cov=matrix(1,1,1,dimnames=list("date", "date")))
mhd3[[1]]
##- Matching within calipers of 3 years
fullmatch(mhd1/(mhd3<3))
}
\keyword{nonparametric }% at least one, from doc/KEYWORDS
