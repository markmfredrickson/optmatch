# Constant to control the maximum feasible (sub)problem
MAX_FEASIBLE <- 1e07

#' (Internal) Sets up the default values for maximum feasible problems
#'
#' @return NULL
setFeasibilityConstants <- function() {
  options("optmatch_warn_on_big_problem" = TRUE)
  options("optmatch_max_problem_size" = MAX_FEASIBLE)
}

#' What is the maximum allowed problem size?
#'
#' To prevent users from starting excessively large matching problems, the
#' maximum problem size is limited by \code{options("optmatch_max_problem_size")}.
#' This function a quick helper to assist fetching this value as a scalar. If the
#' option isn't set, the function falls back to the default value, hard coded in
#' the \code{optmatch} package.
#'
#' @seealso \code{\link{options}}
#' @return logical
#' @examples optmatch:::getMaxProblemSize() > 1 & optmatch:::getMaxProblemSize() < 1e100
#' @export
getMaxProblemSize <- function() {
  tmp <- options("optmatch_max_problem_size")[[1]]
  if (is.null(tmp[[1]])) {
    return(MAX_FEASIBLE)
  }
  return(tmp[[1]])
}

#' Find the minimal exact match factors that will be feasible for a
#' given maximum problem size.
#'
#' The \code{\link{exactMatch}} function creates a smaller matching problem by
#' stratifying observations into smaller groups. For a problem that is larger
#' than maximum allowed size, \code{minExactMatch} provides a way to find the
#' smallest exact matching problem that will allow for matching.
#'
#' \code{x} is a formula of the form \code{Z ~ X1 + X2}, where
#' \code{Z} is indicates treatment or control status, and \code{X1} and \code{X2} are variables
#' can be converted to factors. Any additional arguments are passed to \code{\link{model.frame}}
#' (e.g., a \code{data} argument containing \code{Z}, \code{X1}, and \code{X2}).
#'
#' The the arguments \code{scores} and \code{width} must be passed together.
#' The function will apply the caliper implied by the scores and the width while
#' also adding in blocking factors.
#'
#' @param x The object for dispatching.
#' @param scores Optional vector of scores that will be checked against a caliper width.
#' @param width Optional width of a caliper to place on the scores.
#' @param maxarcs The maximum problem size to attempt to fit.
#' @param ... Additional arguments for methods.
#' @return A factor grouping units, suitable for \code{\link{exactMatch}}.
#' @export
minExactMatch <- function(x, scores = NULL, width = NULL, maxarcs = 1e07, ...) {

  if (length(x) < 3) {
    stop("Formula must be of the form Z ~ X1 + X2 + ...")
  }

  if ((is.null(scores) && !is.null(width) || (!is.null(scores) && is.null(width)))) {
    stop("You must pass both 'scores' and 'width' arguments")
  }

  parts <- rownames(attr(terms(x), "factors"))

  lhs <- parts[1]
  rhs <- parts[-1]
  k <- length(rhs)

  bigzb <- fmla2treatedblocking(x, ...)

  unblockedarcs <- sum(bigzb$Z) * sum(1 - bigzb$Z)
  if (unblockedarcs < maxarcs) {
    return(as.factor(rep(1, dim(bigzb)[1])))
  }

  msg <- getOption("optmatch_verbose_messaging")
  if (msg) {
    warning("minExactMatch: problem is large enough to require blocking. Entering loop.", date())
  }

  previous <- rep(NA, dim(bigzb)[1]) # we store good subgroups here

  for(i in 1:k) {
    fmla <- as.formula(paste(lhs, "~", paste(rhs[1:i], collapse = "+")),
                       env = attr(x, ".Environment"))

    # generate a data.frame of treatment status and the blocking factor:
    z.b <- fmla2treatedblocking(fmla, ...)

    # we create a new factor that includes all previous levels.
    B <- factor(z.b$B, levels = c(levels(previous), levels(z.b$B)))
    B[!is.na(previous)] <- previous[!is.na(previous)]
    B <- factor(B) # toss unused levels

    if (!is.null(scores)) {
      arcs <- caliperSize(scores, z.b$Z, width, structure = B)
    } else {
      arcs <- tapply(z.b$Z, list(B), function(grp) { sum(grp) * sum(1 - grp) })
    }

    good <- arcs <= maxarcs

    if (all(good[!is.na(good)])) { # some levels may be NAs
        if (msg) {
            warning("minExactMatch: exiting loop. Arcs:", arcs, "Selected levels:" , levels(B), date())
        }

      names(B) <- rownames(z.b)
      return(B)
    }

    previous <- factor(B, levels = names(arcs)[good])
  }

  stop("Unable to create sufficiently small problem. Please provide more stratifying variables.")
}

#' (Internal) A helper function to turn formulas into treatment and blocking variables
#'
#' Given a function and any of the arguments normally passed to model.frame,
#' this function will return a data.frame with two columns: a treatment indicator
#' and a blocking factor.
#'
#' @param x A formula
#' @param ... Arguments to be passed to model.frame (e.g. \code{data})
#' @return data.frame containing two columns: \code{Z} is a treatment indicator,
#' \code{B} is a blocking factor
fmla2treatedblocking <- function(x, ...) {

  mf <- model.frame(x, ...)

  blocking <- interaction(mf[,-1])
  treatment <- toZ(mf[,1])

  df <- data.frame(Z = treatment, B = blocking)
  return(df)
}

#' (Internal) Determines how many other units fall within a caliper distance
#'
#' The matching functions \code{\link{fullmatch}} and \code{\link{pairmatch}}
#' have a maximum problem size, based on the number of comparisons between treated
#' and control units. For a completely dense problem, in which every treated units
#' is compared to every control unit there are \code{length(treated) *
#' length(control)} comparisons. A caliper restricts which comparisons are valid,
#' disallowing matches of treated and control pairs that are too far apart. A
#' caliper can significantly decrease the size of a matching problem. The
#' \code{caliperSize} function reports exactly who many valid treated-control
#' comparisons remain after applying a caliper of the given width.
#'
#' @param scores A numeric vector of scores providing 1-D position of units
#' @param z Treatment indicator vector
#' @param width Width of caliper, must be positive
#' @param structure Grouping factor to use in computation
#' @return numeric Total number of pairwise distances remaining after the caliper is placed.
caliperSize <- function(scores, z, width, structure = NULL) {
  if (width <= 0) {
    stop("Invalid caliper width. Width must be positive.")
  }

  if (is.null(structure)) {
    z <- toZ(z)
    treated <- scores[z]
    control <- scores[!z]

    # the following uses findInterval, which requires a sorted vector
    # there may be a speed increase in pulling out the guts of that function and calling them directly
    control <- sort(control)

    stops <- findInterval(treated + width + .Machine$double.eps, control)
    starts <- length(control) - findInterval(-(treated - width -
                                               .Machine$double.eps), rev(-control))
    return(sum(stops - starts))
  }

  # structure is supplied. split up the problem in to blocks and solve those
  results <- sapply(split(data.frame(scores, z), structure), function(x) { caliperSize(x$scores, x$z, width) })
  return(results)

}
# small <- data.frame(y = sample.int(100, 1000, replace = T), z = rep(c(1,0), 500))
# medium <- data.frame(y = sample.int(1000,10000, replace = T), z = rep(c(1,0), 5000))
# big <- data.frame(y = sample.int(10000,100000, replace = T), z = rep(c(1,0), 50000))
#
# system.time(optmatch:::caliperSize(small$y, small$z, 10))
# system.time(optmatch:::caliperSize(medium$y, medium$z, 100))
# system.time(optmatch:::caliperSize(big$y, big$z, 1000))
#
# navie algorithm timings
# 1.6 GHZ Intel Core Duo, 2 GB RAM (used one process, memory did not seem to be an issue)
#
# > system.time(optmatch:::caliperSize(small$y, small$z, 10))
#    user  system elapsed
#   0.030   0.003   0.034
# > system.time(optmatch:::caliperSize(medium$y, medium$z, 100))
#    user  system elapsed
#   2.259   0.011   2.286
# > system.time(optmatch:::caliperSize(big$y, big$z, 1000))
# ^C ... this was taking a long time.
# Timing stopped at: 171.166 29.095 209.98
#
# scaling is much worse than linear. Probably n^2.
#
# New version, 2 core 2GHZ i7, 8GB RAM
# > system.time(optmatch:::caliperSize(small$y, small$z, 10))
#    user  system elapsed
#   0.013   0.002   0.016
# > system.time(optmatch:::caliperSize(medium$y, medium$z, 100))
#    user  system elapsed
#   0.332   0.001   0.333
# > system.time(optmatch:::caliperSize(big$y, big$z, 1000))
#    user  system elapsed
#  29.820   1.584  31.406
#
# Not exactly linear, but not bad!

#' (Internal) Returns a reasonable upper bound on the arcs remaining after placing a caliper.
#'
#' @param scores A numeric vector of scores providing 1-D position of units
#' @param z Treatment indicator vector
#' @param width Width of caliper, must be positive.
#' @param structure Optional factor variable that groups the scores, as would
#' be used by \code{\link{exactMatch}}. Including structure allows for wider
#' calipers.
#' @importFrom graphics hist
#' @return numeric Total number of pairwise distances remaining after the caliper is placed.
caliperUpperBound <- function(scores, z, width, structure = NULL) {

  if (width <= 0) {
    stop("Invalid caliper width. Width must be positive.")
  }

  if (is.null(structure)) {
    z <- toZ(z)
    bins <- seq(min(scores) - width, max(scores) + width, by = width)
    treated <- scores[z]
    control <- scores[!z]

    # the `cut` docs indicate this funciton is faster to get the counts
    control.counts <- hist(control, breaks = bins, plot = FALSE)$counts

    where.treated <- findInterval(treated, bins)
    upper.bound <- 0
    for (i in where.treated) {
      a <- control.counts[i - 1]
      b <- control.counts[i + 1]
      upper.bound <- upper.bound + control.counts[i] + ifelse(is.na(a), 0, a) + ifelse(is.na(b), 0, b)
    }

    return(upper.bound)
  }

  # structure is not null if we get this far
  results <- sapply(split(data.frame(scores, z), structure), function(x) { caliperUpperBound(x$scores, x$z, width) })
  return(sum(results))
}

#' Find the maximum caliper width that will create a feasible problem.
#'
#' Larger calipers permit more possible matches between treated and control
#' groups, which can be better for creating matches with larger effective sample
#' sizes. The downside is that wide calipers may make the matching problem too big
#' for processor or memory constraints. \code{maxCaliper} attempts to find a
#' caliper value, for a given vector of scores and a treatment indicator, that
#' will be possible given the maximum problem size constraints imposed by
#' \code{\link{fullmatch}} and \code{\link{pairmatch}}.
#'
#' @param scores A numeric vector of scores providing 1-D position of units
#' @param z Treatment indicator vector
#' @param widths A vector of caliper widths to try, will be sorted largest to smallest.
#' @param structure Optional factor variable that groups the scores, as would
#' be used by \code{\link{exactMatch}}. Including structure allows for wider
#' calipers.
#' @param exact A logical indicating if the exact problem size should be
#' computed (\code{exact = TRUE}) or if a more computationally efficient upper
#' bound should be used instead (\code{exact = FALSE}). The upper bound may lead
#' to narrower calipers, even if wider calipers would have sufficed using the
#' exact method.
#' @return numeric The value of the largest caliper that creates a feasible
#' problem. If no such caliper exists in \code{widths}, an error will be
#' generated.
#' @export
maxCaliper <- function(scores, z, widths, structure = NULL, exact = TRUE) {
  widths <- sort(widths, decreasing = T)

  f <- caliperUpperBound
  if (exact) { f <- caliperSize}

  for (w in widths) {
    if(all(f(scores, z, w, structure = structure) <= getMaxProblemSize())) {
      return(w)
    }
  }

 stop("Insufficiently small caliper size. Please consider smaller values.")
}
