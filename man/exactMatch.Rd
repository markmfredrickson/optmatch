\name{exactMatch}
\alias{exactMatch}
\title{Generate an exact matching set of subproblems}
\description{
  An exact match is one based on a factor. Within a level, all observations are allowed to be matched.
  An exact match can be combined with another distance matrix to create a set of matching subproblems.
}
\usage{
  exactMatch(x, ...)
  exactMatch(x, data = NULL, subset = NULL, na.action = NULL)
}
\arguments{
  \item{x}{A factor vector or a formula}
  \item{data}{Pptional \code{data.frame} for the formula method.}
  \item{subset}{Optional subset argument as used by \code{model.matrix}.}
  \item{na.action}{Optional \code{NA} handling argument for \code{model.matrix}}
  \item{...}{Arguments for future methods}
}
\details{
  \code{exactMatch} creates a block diagonal matrix of 0s and \code{Inf}s. The pairs
  with 0 entries are within the same level of the factor and legitimate matches.
  \code{\Inf} indicates units in different levels. \code{exactMatch} replaces the
  \code{structure.fmla} argument to several functions in previous versions of Optmatch.
}
\value{
  Object of class \code{DistanceSpecification}, which is suitable to be given
  as \code{distance} argument to \code{\link{fullmatch}} or
  \code{\link{pairmatch}}. The exact match will be only zeros and \code{Inf} values,
  indicating a possible match or no possible match, respectively. It can be added to
  a another distance matrix to create a subclassed matching problem.
}
\author{Mark M. Fredrickson}


\seealso{\code{\link{caliper}}, \code{\link{mdist}}, \code{\link{fullmatch}}, \code{\link{pairmatch}}}

\examples{
data(nuclearplants)

### First generate a standard propensity score
ppty <- glm(pr~.-(pr+cost), family = binomial(), data = nuclearplants)
ppty.distances <- mdist(ppty)

### Only allow matches within the partial turn key plants
pt.em <- exactMatch(pr ~ pt, data = nuclearplants)
as.matrix(pt.em)

### Blunt matches:
match.pt.em <- fullmatch(pt.em)
print(match.pt.em, grouped = T)

### Combine the propensity scores with the subclasses:
match.ppty.em <- fullmatch(ppty.distances + pt.em)
print(match.ppty.em, grouped = T)


}
\keyword{nonparametric }% at least one, from doc/KEYWORDS
