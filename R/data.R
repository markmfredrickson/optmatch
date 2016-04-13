#' The data relate to the construction of 32 light water reactor (LWR)
#' plants constructed in the U.S.A in the late 1960's and early
#' 1970's.  The data was collected with the aim of predicting the cost
#' of construction of further LWR plants.  6 of the power plants had
#' partial turnkey guarantees and it is possible that, for these
#' plants, some manufacturers' subsidies may be hidden in the quoted
#' capital costs.
#'
#' @source The data were obtained from the \code{boot} package, for
#'   which they were in turn taken from Cox and Snell (1981). Although
#'   the data themselves are the same as those in the \code{nuclear}
#'   data frame in the \code{boot} package, the row names of the data
#'   frame have been changed.  (The new row names were selected to
#'   ease certain demonstrations in \code{optmatch}.)
#'
#' This documentation page is also adapted from the \code{boot}
#' package, written by Angelo Canty and ported to R by Brian Ripley.
#'
#' @title Nuclear Power Station Construction Data
#' @format A data frame with 32 rows and 11 columns
#' \itemize{
#'   \item cost: The capital cost of construction in millions of
#'   dollars adjusted to 1976 base.
#'   \item date: The date on which the construction permit was issued.
#'   The data are measured in years since January 1 1990 to the
#'   nearest month.
#'   \item t1: The time between application for and issue of the
#'   construction permit.
#'   \item t2: The time between issue of operating license and
#'   construction permit.
#'   \item cap: The net capacity of the power plant (MWe).
#'   \item pr: A binary variable where \code{1} indicates the prior
#'   existence of a LWR plant at the same site.
#'   \item ne: A binary variable where \code{1} indicates that the
#'   plant was constructed in the north-east region of the U.S.A.
#'   \item ct: A binary variable where \code{1} indicates the use of a
#'   cooling tower in the plant.
#'   \item bw: A binary variable where \code{1} indicates that the
#'   nuclear steam supply system was manufactured by Babcock-Wilcox.
#'   \item cum.n: The cumulative number of power plants constructed by
#'   each architect-engineer.
#'   \item pt: A binary variable where \code{1} indicates those plants
#'   with partial turnkey guarantees.
#' }
#' @references Cox, D.R. and Snell, E.J. (1981) \emph{Applied Statistics:
#'   Principles and Examples}. Chapman and Hall.
#' @keywords datasets
"nuclearplants"

#' This matrix gives discrepancies between light water nuclear power
#' plants of two types, seven built on the site of an existing plant
#' and 19 built on new sites.  The discrepancies summarize differences
#' in two covariates that are predictive of the cost of building a
#' plant.
#'
#' @source The data appear in Cox, D.R. and Snell, E.J. (1981),
#'   \emph{Applied Statistics: Principles and Examples}, p.82 (Chapman
#'   and Hall), and are due to W.E. Mooz.
#'
#' @title Dissimilarities of Some U.S. Nuclear Plants
#' @format A matrix with 7 rows and 19 columns
#' @keywords dataset
#' @references Rosenbaum, P.R. (2002), \emph{Observational Studies},
#'   Second Edition, p.307 (Springer).
"plantdist"
