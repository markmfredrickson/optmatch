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

