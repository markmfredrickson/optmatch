### fill.NAs parses out either a data frame or a formula
### and then hands each column over to the appropriate fill.column
### method. Missingness indicators are added for each column that is
### imputed.

fill.NAs <- function(x, data = NULL, all.covs=FALSE, contrasts.arg=NULL) {
  # if x is present alone, it must be a data.frame
  # if x is present with data, x is a formula
  # all.covs is logical indicating whether or not we want a response variable

  if (is.null(data)) {
    if (!is.data.frame(x) && !is.matrix(x)) {
      stop("Single argument must be a data frame")
    }

    data <- as.data.frame(x) # in case it is a matrix

    # swap the arguments around
    if(!all.covs){
	    response <- colnames(data)[1] # the name of the response var

	    x <- as.formula(paste0("~. - `", response, "`")) # a quick hack to make a formula that is the data, should probably eval col names
      # Added ticks in case response has "+" or other formula characters in it
    }
    if(all.covs){
	    x <- as.formula(~.)
    }

  } else {
    if (inherits(x, "formula")) {
      if(!is.data.frame(data) && !is.matrix(data)) {
        stop("Data argument required for formulas")
      }
      data <- as.data.frame(data) # in case it is a matrix
      ttt <- terms(x, data = data)

      response <- all.vars(ttt)[attr(ttt, "response")]
      data <- data[all.vars(ttt)]
    }
    # TODO should be one more error condition here, if neither formula or d.f
  }

  # not.response <- colnames(data)[colnames(data) != response]

  # find missingness indicators for each column
  original.NAs <- sapply(data, function(i) { any(is.na(i)) })
  original.names <- colnames(data)[original.NAs]


  # create a model matrix from the data; transforms of NA should be NA
  modmat <- model.matrix(x, model.frame(x, data, na.action = na.pass), contrasts.arg=contrasts.arg)
  modmat <- as.data.frame(modmat)
  # remove the intercept, if any
  modmat["(Intercept)"] <- NULL

  # shortcircuit if there are no additional NAs to add
  if(!any(original.NAs)) {
    result <- cbind(data[response], modmat)
    return(result)
  }

  # indicator columns for missing data, only for original missing columns, not transforms
  NA.columns <- sapply(data[original.names], function(column) {
    is.na(column)
  })
  colnames(NA.columns) <- paste(colnames(NA.columns), "NA", sep = ".")

  # of the remaining columns, find those with missingness
  expanded.NAs <- colnames(modmat)[apply(modmat, 2, function(i) { any(is.na(i))})]
  # fill in the columns with missingness
  # NB: fill.column.numeric is hard coded as value of model.matrix is always numeric. no need for a generic fn.
  modmat[expanded.NAs] <- sapply(modmat[expanded.NAs], optmatch:::fill.column.numeric, simplify = F)

  if(!all.covs){
  result <- cbind(data[response], modmat, NA.columns)
  }
  if(all.covs){
	  result <- cbind(modmat,NA.columns)
  }
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
