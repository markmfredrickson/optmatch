% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strata.R
\name{strata}
\alias{strata}
\title{Identify Stratafication Variables}
\usage{
strata(...)
}
\arguments{
\item{...}{any number of variables of the same length.}
}
\value{
the variables with appropriate labels
}
\description{
This is a special function used only in identifying the strata variables
when defining an \code{exactMatch} during a call to \code{fullmatch},
\code{pairmatch}, or \code{match_on}. It should not be called externally.
}
\examples{
data(nuclearplants)
fullmatch(pr ~ cost + strata(pt), data = nuclearplants)
}
