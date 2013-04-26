\name{fill.NAs}
\alias{fill.NAs}
\title{Create missingness indicator variables and non-informatively fill in missing values}
\description{
  Given a \code{data.frame} or \code{formula} and data, \code{fill.NAs()} returns an
  expanded data frame, including a new missingness flag for each variable
  with missing values and replacing each missing entry with a
  value representing a reasonable default for missing values in its column. 
  Functions in the formula are supported, with transformations happening before \code{NA} replacement.
  The expanded data frame is useful for propensity modeling and balance checking when there are
  covariates with missing values. 
}

\usage{
  fill.NAs(x, data = NULL, all.covs = FALSE, contrasts.arg=NULL)

}

\arguments{
  \item{x}{Can be either a data frame (in which case the data argument should be \code{NULL}) or
    a formula (in which case data must be a data.frame)}

  \item{data}{If x is a formula, this must be a data.frame. Otherwise it will be ignored.}

  \item{all.covs}{Should the response variable be imputed? For formula
  \code{x}, this is the variable on the left hand side. For \code{data.frame}
  \code{x}, the response is considered the first column. }

\item{contrasts.arg}{(from \code{model.matrix}) A list, whose entries are values (numeric
    matrices or character strings naming functions) to be used
    as replacement values for the \code{\link{contrasts}}
    replacement function and whose names are the names of
    columns of \code{data} containing \code{\link{factor}}s.}
}

\details{
  \code{fill.NAs} prepares data for use in a model or matching procedure by filling in
  missing values with minimally invasive substitutes. Fill-in is performed column-wise,
  with each column being treated individually. For each column that is missing, a new column
  is created of the form \dQuote{ColumnName.NA} with indicators for each observation that is missing
  a value for \dQuote{ColumnName}.

  The replacement value used to fill in a missing value is simple mean replacement. For transformations
  of variables, e.g. \code{y ~ x1 * x2}, the transformation occurs first. The transformation column will be
  \code{NA} if any of the base columns are \code{NA}. Fill-in occurs next, replacing all missing values
  with the observed column mean. This includes transformation columns.

  Data can be passed to \code{fill.NAs} in two ways. First, you can simply pass a \code{data.frame}
  object and \code{fill.NAs} will fill every column. Alternatively, you can pass a \code{formula} and
  a \code{data.frame}. Fill-in  will only be applied to columns specifically used in the formula. Prior to
  fill-in, any functions in the formula will be expanded. If any arguments to the functions are \code{NA},
  the function value will also be \code{NA} and subject to fill-in.

  By default, \code{fill.NAs} does not impute the response variable. This is to
  encourage more sophisticated imputation schemes when the response is a
  treatment indicator in a matching problem. This behavior can be overridden by
  setting \code{all.covs = TRUE}.
}

\value{
  A \code{data.frame} with all \code{NA} values replaced with mean values
  and additional indicator columns for each column including
  missing values. Suitable for directly passing to \code{\link{lm}} or other
  model building functions to build propensity scores. 
}
\author{ Mark M. Fredrickson and Jake Bowers }

\references{
  von Hipple, Paul T. (2009) \sQuote{How to impute interactions, squares, and other transformed variables,}
    \emph{Sociological Methodlogy}, \bold{39}(1), 265 -- 291.}

\seealso{\code{\link{match_on}}, \code{\link{lm}}}

\examples{

data(nuclearplants)
### Extract some representative covariates:
np.missing <- nuclearplants[c('t1', 't2', 'ne', 'ct', 'cum.n')]

### create some missingness in the covariates
n <- dim(np.missing)[1]
k <- dim(np.missing)[2]

for (i in 1:n) {
  missing <- rbinom(1, prob = .1, size = k)
  if (missing > 0) {
    np.missing[i, sample(k, missing)] <- NA      
  }
}

### Restore outcome and treatment variables:
np.missing <- data.frame(nuclearplants[c('cost', 'pr')], np.missing)

### Fit a propensity score but with missing covariate data flagged 
### and filled in, as in Rosenbaum and Rubin (1984, Appendix):
(np.glm <- glm(fill.NAs(pr ~ t1 * t2, data=np.missing),
family=binomial))

# the previous call is equivalent to:
# glm(pr ~ t1 + t2 + `t1:t2` + t1.NA + t2.NA, fill.NAs(np.missing), family =
#  binomial())

### produce a matrix of propensity distances based on the propensity model 
### with fill-in and flagging. Then perform pair matching on it:
pairmatch(match_on(np.glm))

## fill NAs without using treatment contrasts by making a list of contrasts for each factor
## following hints from http://stackoverflow.com/a/4569239/161808

np.missing$t1F<-factor(np.missing$t1)
cov.factors <- sapply(np.missing[,c("t1F","t2")],is.factor) 
cov.contrasts <- lapply(np.missing[,names(cov.factors)[cov.factors],drop=FALSE],contrasts,contrasts=FALSE)
## make a data frame filling the missing covariate values, but without excluding any levels of any factors
np.noNA2<-fill.NAs(pr~t1F+t2,data=np.missing,contrasts.arg=cov.contrasts)


}

\keyword{manip}% at least one, from doc/KEYWORDS
