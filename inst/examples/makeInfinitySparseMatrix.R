example.matrix <- matrix(c(1,2,Inf, 3,Inf,4, Inf,Inf,Inf), byrow = TRUE, nrow = 3,
                          dimnames = list(letters[1:3], LETTERS[24:26]))

optmatch:::as.InfinitySparseMatrix(example.matrix)

# create the same sparse matrix directly, function will create the appropriate dims
# the data are in a different order, but the indices are correct
(example.ism <- 
  optmatch:::makeInfinitySparseMatrix(c(1,2,3,4), 
                                      c(1,2,1,3), 
                                      c(1,1,2,2), 
                                      LETTERS[24:26],
                                      letters[1:3]))
  

