### fill.NAs parses out either a data frame or a formula
### and then hands each column over to the appropriate fill.column
### method. Missingness indicators are added for each column that is
### imputed.

#' Given a \code{data.frame} or \code{formula} and data,
#' \code{fill.NAs()} returns an expanded data frame, including a new
#' missingness flag for each variable with missing values and
#' replacing each missing entry with a value representing a reasonable
#' default for missing values in its column.  Functions in the formula
#' are supported, with transformations happening before \code{NA}
#' replacement.  The expanded data frame is useful for propensity
#' modeling and balance checking when there are covariates with
#' missing values.
#'
#' \code{fill.NAs} prepares data for use in a model or matching
#' procedure by filling in missing values with minimally invasive
#' substitutes. Fill-in is performed column-wise, with each column
#' being treated individually. For each column that is missing, a new
#' column is created of the form \dQuote{ColumnName.NA} with
#' indicators for each observation that is missing a value for
#' \dQuote{ColumnName}.  Rosenbaum and Rubin (1984, Sec. 2.4 and
#' Appendix B) discuss propensity score models using this data
#' structure.
#'
#' The replacement value used to fill in a missing value is simple
#' mean replacement. For transformations of variables, e.g. \code{y ~
#' x1 * x2}, the transformation occurs first. The transformation
#' column will be \code{NA} if any of the base columns are
#' \code{NA}. Fill-in occurs next, replacing all missing values with
#' the observed column mean. This includes transformation columns.
#'
#' Data can be passed to \code{fill.NAs} in two ways. First, you can
#' simply pass a \code{data.frame} object and \code{fill.NAs} will
#' fill every column. Alternatively, you can pass a \code{formula} and
#' a \code{data.frame}. Fill-in will only be applied to columns
#' specifically used in the formula. Prior to fill-in, any functions
#' in the formula will be expanded. If any arguments to the functions
#' are \code{NA}, the function value will also be \code{NA} and
#' subject to fill-in.
#'
#' By default, \code{fill.NAs} does not impute the response
#' variable. This is to encourage more sophisticated imputation
#' schemes when the response is a treatment indicator in a matching
#' problem. This behavior can be overridden by setting \code{all.covs
#' = TRUE}.
#'
#' @title Create missingness indicator variables and non-informatively fill in missing values
#'
#' @param x Can be either a data frame (in which case the data
#'   argument should be \code{NULL}) or a formula (in which case data
#'   must be a data.frame)
#' @param data If x is a formula, this must be a data.frame. Otherwise
#'   it will be ignored.
#' @param all.covs Should the response variable be imputed? For
#'   formula \code{x}, this is the variable on the left hand side. For
#'   \code{data.frame} \code{x}, the response is considered the first
#'   column.
#' @param contrasts.arg (from \code{model.matrix}) A list, whose
#'   entries are values (numeric matrices or character strings naming
#'   functions) to be used as replacement values for the
#'   \code{\link{contrasts}} replacement function and whose names are
#'   the names of columns of \code{data} containing
#'   \code{\link{factor}}s.
#'
#' @return A \code{data.frame} with all \code{NA} values replaced with
#'   mean values and additional indicator columns for each column
#'   including missing values. Suitable for directly passing to
#'   \code{\link{lm}} or other model building functions to build
#'   propensity scores.
#'
#' @author Mark M. Fredrickson and Jake Bowers
#'
#' @references
#' Rosenbaum, Paul R. and Rubin, Donald B. (1984) \sQuote{Reducing
#'   Bias in Observational Studies using Subclassification on the
#'   Propensity Score,} \emph{Journal of the American Statistical
#'   Association}, \bold{79}, 516 -- 524.
#'
#' Von Hipple, Paul T. (2009) \sQuote{How to impute interactions,
#'   squares, and other transformed variables,} \emph{Sociological
#'   Methodology}, \bold{39}(1), 265 -- 291.
#'
#' @seealso \code{\link{match_on}}, \code{\link{lm}}
#'
#' @examples
#' data(nuclearplants)
#' ### Extract some representative covariates:
#' np.missing <- nuclearplants[c('t1', 't2', 'ne', 'ct', 'cum.n')]
#'
#'  ### create some missingness in the covariates
#'  n <- dim(np.missing)[1]
#'  k <- dim(np.missing)[2]
#'
#'  for (i in 1:n) {
#'    missing <- rbinom(1, prob = .1, size = k)
#'    if (missing > 0) {
#'      np.missing[i, sample(k, missing)] <- NA
#'    }
#'  }
#'
#' ### Restore outcome and treatment variables:
#' np.missing <- data.frame(nuclearplants[c('cost', 'pr')], np.missing)
#'
#' ### Fit a propensity score but with missing covariate data flagged
#' ### and filled in, as in Rosenbaum and Rubin (1984, Appendix):
#' np.filled <- fill.NAs(pr ~ t1 * t2, np.missing)
#' # Look at np.filled to establish what missingness flags were created
#' head(np.filled)
#' (np.glm <- glm(pr ~ ., family=binomial, data=np.filled))
#' (glm(pr ~ t1 + t2 + `t1:t2` + t1.NA + t2.NA,
#'                 family=binomial, data=np.filled))
#' # In a non-interactive session, the following may help, as long as
#' # the formula passed to `fill.NAs` (plus any missingness flags) is
#' # the desired formula for the glm.
#' (glm(formula(terms(np.filled)), family=binomial, data=np.filled))
#'
#' ### produce a matrix of propensity distances based on the propensity model
#' ### with fill-in and flagging. Then perform pair matching on it:
#' pairmatch(match_on(np.glm, data=np.filled), data=np.filled)
#'
#' ## fill NAs without using treatment contrasts by making a list of contrasts for
#' ## each factor ## following hints from http://stackoverflow.com/a/4569239/161808
#'
#' np.missing$t1F<-factor(np.missing$t1)
#' cov.factors <- sapply(np.missing[,c("t1F","t2")],is.factor)
#' cov.contrasts <- lapply(
#'   np.missing[,names(cov.factors)[cov.factors],drop=FALSE],
#'   contrasts, contrasts = FALSE)
#'
#' ## make a data frame filling the missing covariate values, but without
#' ## excluding any levels of any factors
#' np.noNA2<-fill.NAs(pr~t1F+t2,data=np.missing,contrasts.arg=cov.contrasts)
#'
#' @keywords manip
#'
#' @export
fill.NAs <- function(x, data = NULL, all.covs = FALSE, contrasts.arg=NULL) {
  # if x is present alone, it must be a data.frame
  # if x is present with data, x is a formula
  # all.covs is logical indicating whether or not we want a response variable

  if (is.null(data)) {
    if (!is.data.frame(x) && !is.matrix(x)) {
      stop("Single argument must be a data frame")
    }

    # swap the arguments around
    data     <- as.data.frame(x) # in case it is a matrix
	  response <- colnames(data)[1] # the name of the response var

    if(!all.covs && dim(data)[2] > 1){
      # Added ticks in case response has "+" or other formula characters in it
	    x <- terms(as.formula(paste0("`", response, "` ~  .")), data = data)
    } else {
	    x <- terms(as.formula(~ .), data = data) # everything, including the response should be imputed
    }
  } else {
    if (inherits(x, "formula")) {
      if(!is.data.frame(data) && !is.matrix(data)) {
        stop("Data argument required for formulas")
      }
    }
  }

  withStrata <- findStrata(x, data)

  data <- as.data.frame(data) # in case it is a matrix
  ttt <- terms(withStrata$newx, data = data)

  response <- all.vars(ttt)[attr(ttt, "response")]
  data.trimmed <- data[all.vars(ttt)]

  # not.response <- colnames(data)[colnames(data) != response]

  # find missingness indicators for each column
  original.NAs <- sapply(data.trimmed, function(i) { any(is.na(i)) })
  original.names <- colnames(data.trimmed)[original.NAs]

  # create a model matrix from the data.trimmed; transforms of NA should be NA
  modmat <- model.matrix(withStrata$newx, model.frame(withStrata$newx, data.trimmed, na.action = na.pass), contrasts.arg=contrasts.arg)
  modmat <- as.data.frame(modmat)
  # remove the intercept, if any
  modmat["(Intercept)"] <- NULL

  # shortcircuit if there are no additional NAs to add
  # if(!any(original.NAs)) {
  #   result <- cbind(data.trimmed[response], modmat)
  #   return(result)
  # }
  result <- modmat

  if(any(original.NAs)) {

    # indicator columns for missing data.trimmed, only for original missing columns, not transforms
    NA.columns <- sapply(data.trimmed[original.names], function(column) {
      is.na(column)
    })
    colnames(NA.columns) <- paste(colnames(NA.columns), "NA", sep = ".")

    # of the remaining columns, find those with missingness
    expanded.NAs <- colnames(modmat)[apply(modmat, 2, function(i) { any(is.na(i))})]
    # fill in the columns with missingness
    # NB: fill.column.numeric is hard coded as value of model.matrix is always numeric. no need for a generic fn.
    if (length(withStrata$strata) > 0 ) {
      sformula <- as.formula(paste("~", paste(withStrata$strata, collapse = "+")))
      tmp <- interaction(model.frame(sformula, data = data, na.action = na.pass))
      for (l in levels(tmp)) {
        idx <- tmp == l & !is.na(tmp)
        modmat[idx, expanded.NAs] <- sapply(modmat[expanded.NAs][idx, , drop = FALSE], fill.column.numeric, simplify = F)
      }
    } else {
      modmat[expanded.NAs] <- sapply(modmat[expanded.NAs], fill.column.numeric, simplify = F)
    }
    result <- cbind(modmat, NA.columns)
  }

  if(!all.covs){
    result <- cbind(data.trimmed[response], result)
    newfmla <- formula(result)
  } else {
    p <- paste0("`", names(result), "`")
    p <- gsub("``", "`", p)
    newfmla <- formula(paste("~", paste(p, collapse="+")))
  }

  if (length(withStrata$strata) > 0) {
    sformula      <- as.formula(paste("~", paste(withStrata$strata, collapse = "+")))
    tmp           <- model.frame(sformula, data = data, na.action = na.pass)
    colnames(tmp) <- all.vars(sformula)
    result        <- cbind(result, tmp)
    newfmla <- update(newfmla, paste0(".~.+", paste(withStrata$strata, collapse = "+")))
  }


  attr(result, "terms") <- terms(newfmla, data = data, specials = "strata")

  return(result)
}


### Column imputation: takes a column and returns filled in values
### Adding NA indicator columns happens elsewhere
### NB: the only that will probably ever be called in numeric now
### that model.matrix gets called ahead of column filling

fill.column <- function(column) {
  UseMethod("fill.column", column)
}

fill.column.numeric <- function(column) {
  nas <- is.na(column)
  cm <- mean(column[!nas])
  column[nas] <- cm
  return(column)
}

fill.column.logical <- function(column) {
  nas <- is.na(column)
  cm <- mean(column[!nas]) > .5
  column[nas] <- cm
  return(column)

}

fill.column.factor <- function(column) {
  # following RITools' imputation function, this adds a level
  # another option would be imputing based on the modal factor
  levels(column) <- c(levels(column),'.NA')
  column[is.na(column)] <- ".NA"
  return(column)
}

fill.column.ordered <- function(column) {

# uses a median imputation scheme -- finds the middle most
# level and uses that
  sorted <- sort(column, NA.last = NA)
  imputed.value <- sorted[floor(length(sorted))/2]
  column[is.na(column)] <- imputed.value
  return(column)

}
