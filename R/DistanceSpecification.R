################################################################################
### Objects for specifying distances, methods to manipulate them
################################################################################

# The following objects are valid DistanceSpecifications, that is, they layout
# how treatment and control units are related.

setClassUnion("DistanceSpecification", c("matrix", "matrix.csr"))

### prepareMatching: DistanceSpecification -> arcs
### where arcs is a data.frame with 3 columns: control, treatment, distance

setGeneric("prepareMatching", function(distances)
  standardGeneric("prepareMatching"))

setMethod("prepareMatching", "matrix",
function(distances) {
  # note: treatment units are the rows, control units are the columns
  cs <- colnames(distances)
  rs <- rownames(distances)
  raw.distances <- as.vector(distances)
  
  control.reps <- length(raw.distances) / length(cs)
  treatment.reps <- length(raw.distances) / length(rs)

  controls <- rep(cs, each = control.reps)
  treatments <- rep(rs, treatment.reps)

  tmp <- data.frame(control = controls, treatment = treatments, distance =
    raw.distances)
  tmp <- tmp[is.finite(tmp$distance),]

  return(tmp) 
})

setMethod("prepareMatching", "matrix.csr",
function(distances) {
  # From the SparseM vignette on the internal structure of the matrix.csr
  # class:

  ### ra: a real array of nnz elements containing the non-zero elements of A,
  ### stored in row order.  Thus, if i < j, all elements of row i precede elements
  ### from row j. The order of elements within the rows is immaterial.

  ### ja: an integer array of nnz elements containing the column indices of the
  ### elements stored in ra.  
  
  ### ia: an integer array of n + 1 elements containing pointers to the beginning
  ### of each row in the arrays ra and ja. Thus ia[i] indicates the position in
  ### the arrays ra and ja where the ith row begins. The last (n + 1)st element of
  ### ia indicates where the n + 1 row would start, if it existed.

  # To get to the cannonical representation, ra is the distances, ja the
  # control column. Treatment requires a little rejiggering.

  rowids <- 1:(distances@dimension[1]) # 1, 2, ... max row number
  treatments <- rep(rowids, diff(distances@ia))
  tmp <- data.frame(control = distances@ja, treatment = treatments, distance =
    distances@ra)
  return(tmp) 
})
