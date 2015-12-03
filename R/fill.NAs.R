### fill.NAs parses out either a data frame or a formula
### and then hands each column over to the appropriate fill.column
### method. Missingness indicators are added for each column that is
### imputed.

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
  newfmla <- x
  
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

    newfmla <- update(newfmla, as.formula(paste0("~ . + ", paste0("`", colnames(NA.columns), "`", collapse = "+"))))
  }

  if(!all.covs){
    result <- cbind(data.trimmed[response], result)
  }

  if (length(withStrata$strata) > 0) {
    sformula      <- as.formula(paste("~", paste(withStrata$strata, collapse = "+")))
    tmp           <- model.frame(sformula, data = data, na.action = na.pass) 
    colnames(tmp) <- all.vars(sformula)
    result        <- cbind(result, tmp)
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
