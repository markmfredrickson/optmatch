################################################################################
# exactMatch: create InfinitySparseMatrices from factors
################################################################################

exactMatch <- function(split, treatment) {

  if (length(split) != length(treatment)) {
    stop("Splitting vector and treatment vector must be same length")  
  }
  
  # defensive programming
  split <- as.factor(split)
  treatment <- as.factor(treatment)

  if (length(levels(treatment)) != 2) {
    stop("Treatment must be a factor of two levels")  
  }

  # the upper level is the treatment condition
  splitT <- split[treatment == 1]
  splitC <- split[!treatment == 0]

  csForTs <- lapply(splitT, function(t) {
    which(t == splitC)
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
}
