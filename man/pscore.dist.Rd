\name{pscore.dist}
\alias{pscore.dist}
\title{Assemble propensity distances and prepare for matching}
\description{
  Extracts scores from a fitted propensity scoring model, 
  assembling them into a discrepancy matrix (or matrices) from which
  \code{pairmatch()} or \code{fullmatch()} can determine optimal matches.
}
\usage{
pscore.dist(glmobject, structure.fmla = NULL, standardization.scale=sd)
}
\arguments{
  \item{glmobject}{ A fitted propensity modeling object, produced by a
    call to \code{glm()} or, say, \code{bayesglm} (from package
    \dQuote{arm}) or \code{brglm()} (from the \dQuote{brglm} package). }
  \item{structure.fmla}{ Optional formula argument specifying subclasses
  within which matches are to be performed.  If omitted, no
  subclassification is done. If it is given, the RHS of this formula gives
  variables on which to stratify the sample prior to matching. }
  \item{standardization.scale}{Scalar-valued function or \code{NULL},
    defaulting to \code{sd}.}
}
\details{
  \code{glmobject} need not necessarily be the result of a call to
  \code{glm}.  It should be a list with elements: \code{y}, a vector that
  is positive for treatment subjects and nonpositive for controls;
  \code{linear.predictors}, containing the propensity scores; and
  \code{data}, the data frame from which propensities were made.
  
  The purpose of giving a \code{structure.fmla} argument is to speed
  up large problems.  Variables appearing on its
  right-hand side will be interacted to create these subclasses; the
  same variables should also have appeared on the RHS of the formula
  used to specify the propensity model.

  If non-null, argument
    \code{standardization.scale} should be a scalar-valued function of a
    vector argument.  It is applied separately to
    treatment and control propensity scores to determine their scale;
    propensity distances will be scaled by a pooling of these two values.
    If \code{NULL}, no scaling is supplied, and distances are in terms
    of the linear propensity score, that is the logits of the
    conditional probabilities. (Another good choice is \code{\link{mad}},
    which is robust to outliers but gives results similar to \code{sd}
    in the absence of outliers.)
    
  Underneath, the function makes use of \code{\link{makedist}}. So it keeps track of
  metadata useful for matching, as described on the help page for that function.
}
\value{
  Object of class \code{optmatch.dlist}, which is suitable to be given
  as \code{distance} argument to \code{\link{fullmatch}} or
  \code{\link{pairmatch}}.

  Specifically, a list of matrices, one for each subclass defined by the
  interaction of variables appearing on the right hand side of
  \code{structure.fmla}.  Each of these is a number of treatments by
  number of controls matrix of propensity distances.  The distances are
  differences of the linear predictor from the propensity model, rather than
  differences of estimated probabilities, avoiding compression of
  estimated propensities near 0 and 1 (Rosenbaum and Rubin, 1985).   They will
  have been scaled by the pooled SD of propensity scores in the
  treatment and control groups, so that a caliper of .25 pooled SDs on
  the propensity score can be coded as \code{value/(value<=.25)}; see
  the examples.  
  The list also carries some metadata as attributes, data that are not of direct interest to
  the user but are useful to  \code{fullmatch()} and \code{pairmatch()}.
}
\author{ Ben B. Hansen}

\references{
  P.~R. Rosenbaum and D.~B. Rubin (1985),
\sQuote{Constructing a control group using multivariate matched sampling
  methods that incorporate the propensity score},
 \emph{The American Statistician}, \bold{39} 33--38.
}

\seealso{ \code{\link{makedist}}, \code{\link{mahal.dist}}, \code{\link{fullmatch}}, \code{\link{pairmatch}},
          \code{\link{caliper}}}
\examples{
data(nuclearplants)
psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
psd <- pscore.dist(psm)
fullmatch(psd)
fullmatch(caliper(.25, psd)) # Propensity matching with calipers 
}
\keyword{nonparametric }% at least one, from doc/KEYWORDS
