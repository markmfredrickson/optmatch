caliper <- function(width, ..., exclude = c(), penalty = Inf) {

  start <- mdist(...) # returns an optmatch.dlist object

  # pseudo:
  # set res widths of same name to start values
  res <- sapply(start, function(mat) {
    mat[mat > width] <- penalty * mat[mat > width]
    mat[!(mat > width)] <- 0
    return(mat)
  }, simplify = F)
  f <- revertrows(exclude)
  res <- mapply(f, res, start, SIMPLIFY = F)

  # after mapplying, res becomes a plain list, make it a dlist
  class(res) <- list("optmatch.dlist", "list")
  attributes(res) <- attributes(start)

  return(res)
}

revertrows <- function(exclude) { 
  function(updateme, orig) {
    rns <- exclude[exclude %in% rownames(orig)]
    cns <- exclude[exclude %in% colnames(orig)]
    updateme[rns, ] <- 0 
    updateme[,cns] <- 0
    return(updateme)
  }
}

# below was a first cut at allowing individual widths per observation
# I scrapped it for the simpler revertrows for now.
#createWidths <- function(width, mat, namedwdiths) {
#  # width is the default width for any non-named widths
#  # mat is a single entry in a optmatch.dlist object
#  # namedwidths is a vector with names 
#  mrows <- dim(mat)[1]
#  mcols <- dim(mat){2}
#  newmat <- matrix(width, nrows = mrows, mcols)
#  for (i in 1:mrows) {
#    for (j in 1:mcols)
#      newmat[i, j] <- min()
#    }
#  }
    
#}
