caliper <- function(width, ..., exclude = c(), compare = `<=`) {

  start <- as.InfinitySparseMatrix(mdist(...)) # returns a DistanceSpecification

  excluded.rows <- which(start@rownames %in% exclude)
  excluded.cols <- which(start@colnames %in% exclude)

  x <- subset(start, compare(start, width) | 
                     start@rows %in% excluded.rows |
                     start@cols %in% excluded.cols)

  x@.Data <- rep(0, length(x@.Data))

  return(x)
  
  # this might be useful if I reinstate the penalty argument
  #
  # res <- sapply(start, function(mat) {
  #   mat[mat > width] <- penalty * mat[mat > width]
  #   mat[!(mat > width)] <- 0
  #   return(mat)
  # }, simplify = F)
  # f <- revertrows(exclude)
  # res <- mapply(f, res, start, SIMPLIFY = F)

  # # after mapplying, res becomes a plain list, make it a dlist
  # class(res) <- list("optmatch.dlist", "list")
  # attributes(res) <- attributes(start)

  # return(res)
}

# revertrows <- function(exclude) { 
#   function(updateme, orig) {
#     rns <- exclude[exclude %in% rownames(orig)]
#     cns <- exclude[exclude %in% colnames(orig)]
#     updateme[rns, ] <- 0 
#     updateme[,cns] <- 0
#     return(updateme)
#   }
# }

