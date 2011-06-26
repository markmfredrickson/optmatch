################################################################################
# exactMatch: create InfinitySparseMatrices from factors
################################################################################

setGeneric("exactMatch", 
  def = function(x, ...) standardGeneric("exactMatch"))

setMethod(exactMatch, "vector", function(x, treatment) {
  if (length(x) != length(treatment)) {
    stop("Splitting vector and treatment vector must be same length")  
  }
  
  # defensive programming
  x <- as.factor(x)
  treatment <- as.factor(treatment)

  if (length(levels(treatment)) != 2) {
    stop("Treatment must be a factor of two levels")  
  }

  # the upper level is the treatment condition
  xT <- x[treatment == 1]
  xC <- x[!treatment == 0]

  csForTs <- lapply(xT, function(t) {
    which(t == xC)
  })

  cols <- unlist(csForTs)
  tmp <- sapply(csForTs, length)
  rows <- rep(1:(length(csForTs)), tmp)

  if (!is.null(names(treatment))) {
    rns <- names(treatment)[treatment == 1]
    cns <- names(treatment)[treatment == 0]
  } else {
    rns <- NULL
    cns <- NULL
  }

  return(makeInfinitySparseMatrix(rep(0, length(rows)), cols = cols, rows =
    rows, rownames = rns, colnames = cns))  
})

setMethod(exactMatch, "formula", function(x, data = NULL) {
  mf <- model.frame(x)
  
  # formula is expected to be Z ~ B, where b is the blocking factor
  # and Z is treatment, Z ~ B1 + B2 ... is also allowed
  exactMatch(interaction(mf[,-1]), mf[,1]) # use the factor based method
})
