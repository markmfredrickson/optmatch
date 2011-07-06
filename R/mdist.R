################################################################################
# Mdist: distance matrix creation functions
################################################################################

setGeneric("mdist", def = function(x, exclusions = NULL, ...)  standardGeneric("mdist"))

setMethod("mdist", "function", function(x, exclusions = NULL, z = NULL, data = NULL, ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }

  theFun <- match.fun(x)

  makedist(z, data, theFun, exclusions)
})


# mdist method: formula
setMethod("mdist", "formula", function(x, exclusions = NULL, data = NULL, subset = NULL, 
                                       s.matrix = NULL, ...) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf) == "x"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  # TODO: check x for correct format Z ~ ...

  if (!is.null(s.matrix)) {
    # TODO: error check the inv.cov matrix to make sure it is safe
    # should match dimension of mf
  } else {
    # default s.matrix is the inverse covariance matrix
    s.matrix <- solve(cov(mf[,-1])) # don't need Z in the cov matrix
  }

  f <- function(treated, control) {
    n <- dim(treated)[1]
    tmp <- numeric(n) 
    for (i in 1:n) {
      tmp[i] <- sqrt(as.matrix(treated[i,] - control[i,]) %*% s.matrix %*% t(as.matrix(treated[i,] - control[i,]))) 
    }
    return(tmp)
  }

  makedist(mf[,1], mf[,-1], f, exclusions)

})

# mdist method: glm
setMethod("mdist", "glm", function(x, exclusions = NULL, standardization.scale = mad, ...)
{
  stopifnot(all(c('y', 'linear.predictors','data') %in% names(x)))
  z <- x$y > 0
  pooled.sd <- if (is.null(standardization.scale)) {
    1 
  } else {
    szn.scale(x$linear.predictors, z ,standardization.scale)
  }

  lp.adj <- x$linear.predictors/pooled.sd

  f <- function(t, c) { abs(t - c) }
  
  makedist(z, lp.adj, f, exclusions)
})

szn.scale <- function(x, Tx, standardizer = mad, ...) {
  sqrt(((sum(!Tx) - 1) * standardizer(x[!Tx])^2 + 
        (sum(!!Tx) - 1) * standardizer(x[!!Tx])^2) / (length(x) - 2))
}

# mdist method: bigglm
setMethod("mdist", "bigglm", function(x, exclusions = NULL, data = NULL, standardization.scale = mad, ...)
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
})


### mdist method: numeric.
### (mdist can't work with numeric vectors at present,
### but it can return an informative error message).

setMethod("mdist", "numeric", function(x, exclusions = NULL, ...)
{

  stop("No mdist method for numerics.
  Consider using mdist(z ~ ps | strata, data = your.data)
  where ps is your numeric vector, z is your treatment assignment,
  and strata (optional) indicates a stratification variable, all
  columns in your.data")

})

# mdist methods for DistanceSpecifications
# apparently the class union is less important than the true
# type, so the numeric method above gets in the way
setMethod("mdist", "InfinitySparseMatrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

setMethod("mdist", "matrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

