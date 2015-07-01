setClass("ImputedFormula",
         contains = "formula",
         representation = list(imputeLHS = "logical"))

setOldClass("ImputedFormula", c("ImputedFormula", "formula"))
# setAs("ImputedFormula", "formula", function(from) from@.Data) ## dummy coercion here

#' Create a special formula that will impute NAs in predictor columns.
#'
#' @param formula A formula object.
#' @param imputeLHS A flag indicating if we should impute the left hand side of the formula (the repsonse). By default, the response is not imputed.
#' @return A special formula object that will properly impute in models or other expansion.
#' @export
fill.NAs <- function(formula, imputeLHS = FALSE) {
  formula <- as(formula, Class = "ImputedFormula")
  formula@imputeLHS = imputeLHS
  return(formula)
}

#' @S3method model.frame ImputedFormula
model.frame.ImputedFormula <- function(formula, na.action = na.pass, ...) {
  # TODO: should strike missing outcome units unless imputeLHS is true
  tmp <- as(formula, "formula")
  mf <- model.frame(tmp, na.action = na.action, ...)
  class(attr(mf, "terms")) <- c("ImputedFormulaTerms", "terms", "formula")
  return(mf)
}

#' @S3method model.matrix ImputedFormulaTerms
model.matrix.ImputedFormulaTerms <- function(object, data, ...) {
  model.matrix.ImputedFormula(object, data, ...)
}

#' @S3method model.matrix ImputedFormula
model.matrix.ImputedFormula <- function(object, data, ...) { 

  if (missing(data) | is.null(data)) {
    stop("Have not implemented fetching data from environment yet")
  }

  withStrata <- findStrata(object, data)
 
  mm <- model.matrix.default(object, data = data, ...)

  original.NAs <- colnames(data)[sapply(data, function(i) { any(is.na(i)) })]

  if(length(original.NAs) > 0) {

    # indicator columns for missing data.trimmed, only for original missing columns, not transforms
    NA.columns <- sapply(data[original.NAs], function(column) {
      as.numeric(is.na(column))
    })
    colnames(NA.columns) <- paste(colnames(NA.columns), "NA", sep = ".")

    # of the remaining columns, find those with missingness
    mm.NAs <- colnames(mm)[apply(mm, 2, function(i) { any(is.na(i))})]

    # fill in the columns with missingness
    # NB: fill.column.numeric is hard coded as value of model.matrix is always numeric. no need for a generic fn.
    # if (length(withStrata$strata) > 0 ) {
    #   sformula <- as.formula(paste("~", paste(withStrata$strata, collapse = "+")))
    #   tmp <- interaction(model.frame(sformula, data = data, na.action = na.pass))
    #   for (l in levels(tmp)) {
    #     idx <- tmp == l & !is.na(tmp)
    #     modmat[idx, expanded.NAs] <- sapply(modmat[expanded.NAs][idx, , drop = FALSE], fill.column.numeric, simplify = F)
    #   }
    # } else {
    #   modmat[expanded.NAs] <- sapply(modmat[expanded.NAs], fill.column.numeric, simplify = F)
    # }
    mm <- cbind(apply(mm, 2, fill.column.numeric), NA.columns)
  }

  return(mm)
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
