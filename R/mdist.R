mdist <- function(x, structure.fmla = NULL, ...) {
  UseMethod("mdist", x)
}

# mdist method: optmatch.dlist
mdist.optmatch.dlist <- function(x, structure.fmla = NULL, ...) {
  return(x)
} # just return the argument

# mdist method: function
# for the function method, both data and structure fmla are required
# api change from makedist: I think it would make more sense to have
# the function take two data.frames: the treatments in this stratum
# and the controls in the stratum. It could then return the matrix of
# mdists, which the rest of the function would markup with rownames
# etc.
mdist.function <- function(x, structure.fmla = NULL, data = NULL, ...) {

  if (is.null(data) || is.null(structure.fmla)) {
    stop("Both data and the structure formula are required for
    computing mdists from functions.")
  }

  theFun <- match.fun(x)
  parsedFmla <- parseFmla(structure.fmla)

  if(is.null(parsedFmla[[1]])) {
    stop("Treatment variable required")
  }

  if((identical(parsedFmla[[2]], 1) && is.null(parsedFmla[3])) ||
     (length(parsedFmla[[2]]) > 1)) {
    stop("Please specify the grouping as either: z ~ grp or z ~ 1 | grp")
  }

  treatmentvar <- parsedFmla[[1]]

  if(is.null(parsedFmla[3][[1]])) { # I swear subscripting is incomprehensible!
    strata <- parsedFmla[[2]]
  } else {
    strata <- parsedFmla[[3]]
  }


  # split up the dataset by parts
  # call optmatch.dlist maker function on parts, names

    # create a function to produce one distance matrix
  doit <- function(data,...) {
     # indicies are created per chunk to split out treatment and controls
     indices <- data[as.character(treatmentvar)] == 1
     treatments <- data[indices,]
     controls <- data[!indices,]
     distances <- theFun(treatments, controls, ...)

     colnames(distances) <- rownames(controls)
     rownames(distances) <- rownames(treatments)

     return(distances)
  }

  if (!(identical(strata, 1))) {
    if(is.factor(eval(strata, data))){
      ss <- eval(strata, data) ##to preserve existing labels/levels
    } else {
      ss <- factor(eval(strata, data), labels = 'm')
    }
    ans <- lapply(split(data, ss), doit,...)
  } else {
    ans <- list(m1 = doit(data,...))
  }

  attr(ans, 'row.names') <- row.names(data)
  class(ans) <- c('optmatch.dlist', 'list')
  return(ans)
}


# mdist method: formula
mdist.formula <- function(x, structure.fmla = NULL, data = NULL, subset=NULL,...) {
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=="x"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")

  if (length(x)==2)
      if (!is.null(structure.fmla) && length(structure.fmla)==3)
          x <- update.formula(structure.fmla, x) else
  stop("If you give mdist() a formula, it needs to have a left hand side\nin order for mdist.formula() to figure out who is to be matched to whom.",
       call.=FALSE)

  if (isThereAPipe(x)) {
      if (!is.null(structure.fmla))
          warning("I see a pipe, a '|', in the formula you gave as primary argument to mdist().\n So I'm using what's to the right of it in lieu of the RHS of the\n'structure.fmla' argument that you've also given.", call.=FALSE)
      parsed <- parseFmla(x)
    # this block occurs if the grouping factor is present
    # e.g. z ~ x1 + x2 | grp
    if (!is.null(unlist(parsed[3]))) {
      xenv <- environment(x)
      x <- as.formula(paste(as.character(parsed[1:2]), collapse = "~"))
      environment(x) <- xenv
      structure.fmla <- as.formula(paste("~", parsed[[3]]))
    }
  }
  mf$formula <-  makeJoinedFmla(x, structure.fmla)
  mf <- eval(mf, parent.frame())

###  return(mf)
  mahal.dist(x, data = mf, structure.fmla = structure.fmla, ...)
}

isThereAPipe <- function(fmla)
{
    inherits(fmla, "formula") && length(fmla) == 3 && length(fmla[[3]])==3 && fmla[[3]][[1]] == as.name("|")
}
# One big formula from a formula plus astructure formula -- for use w/ model.frame
makeJoinedFmla <- function(fmla, structure.fmla)
{
    if (is.null(structure.fmla)) return(fmla)
    stopifnot(inherits(fmla, "formula"), inherits(structure.fmla, "formula"),
              length(fmla)==3)

l <- length(structure.fmla)
structure.fmla[[l]] <- as.call(c(as.name("+"), as.name("."), structure.fmla[[l]]))

update.formula(fmla, structure.fmla)
}

# mdist method: glm
mdist.glm <- function(x, structure.fmla = NULL, standardization.scale=mad, ...)
{
  pscore.dist(x,  structure.fmla = structure.fmla, standardization.scale=standardization.scale, ...)
}

# parsing formulas for creating mdists
parseFmla <- function(fmla) {

  treatment <- fmla[[2]]
  rhs <- fmla[[3]]
  if (length(rhs) == 3 && rhs[[1]] == as.name("|")) {
    covar <- rhs[[2]]
    group <- rhs[[3]]

  } else {
    covar <- rhs
    group <- NULL
  }

  return(c(treatment, covar, group))

}


# mdist method: bigglm
mdist.bigglm <- function(x, structure.fmla = NULL, data = NULL, standardization.scale=mad, ...)
{
  if (is.null(data))
    stop("data argument is required for computing mdists from bigglms")

  if (!is.data.frame(data))
    stop("mdist doesn't understand data arguments that aren't data frames.")

  if (is.null(structure.fmla) | !inherits(structure.fmla, 'formula'))
    stop("structure.fmla argument required with bigglms.
(Use form 'structure.fmla=<treatment.variable> ~ 1'
 for no stratification before matching)")

theps <- predict(x, data, type='link', se.fit=FALSE)
if (length(theps)!=dim(data)[1])
stop("predict.bigglm() returns a vector of the wrong length;
are there missing values in data?")


Data <-  model.frame(structure.fmla, data=data)
treatmentvar <- as.character(structure.fmla[[2]])
pooled.sd <- if (is.null(standardization.scale)) 1 else
szn.scale(theps, Data[[treatmentvar]], standardizer=standardization.scale,...)

Data$tHePs <- theps/pooled.sd

psdiffs <- function(treatments, controls) {
abs(outer(as.vector(treatments$tHePs),
as.vector(controls$tHePs), `-`))
}

mdist(psdiffs, structure.fmla=structure.fmla,
      data=Data)
}


### mdist method: numeric.
### (mdist can't work with numeric vectors at present,
### but it can return an informative error message).

mdist.numeric <- function(x, structure.fmla = NULL, trtgrp=NULL, ...)
{

  stop("No mdist method for numerics.
  Consider using mdist(z ~ ps | strata, data = your.data)
  where ps is your numeric vector, z is your treatment assignment,
  and strata (optional) indicates a stratification variable, all
  columns in your.data")

}
