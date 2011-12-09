\name{makeInfinitySparseMatrix}
\alias{makeInfinitySparseMatrix}
\alias{as.InfinitySparseMatrix}
\title{Creating sparse matching problems}
\description{
  Create \code{InfinitySparseMatrix} distance specifcations. Finite entries indicate possilbe matches,
while infinite entries indicated non-allowed matches. This data type can be more space efficient for sparse
matching problems.
}
\usage{
makeInfinitySparseMatrix(data, cols, rows, colnames = NULL, rownames = NULL, dimension = NULL)
as.InfinitySparseMatrix(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A vector of distances for the finite (allowed) treatment-control pairs.}
  \item{cols}{A vector indicating the column number for each entry in \code{data}.}
  \item{rows}{A vector indicating the row number for each entry in \code{data}.}
  \item{colnames}{A optional character vector with the columns names of the matrix.}
  \item{rownames}{A optional character vector with the row names of the matrix.}
  \item{dimension}{An optional vector giving the dimensions of the matrix, which can be useful for indicating
  matrices with entirely \code{Inf} rows or columns. If supplied with row and columns names, it must match.}
  \item{x}{A matrix to be converted to an \code{InfinitySparseMatrix}.}
}
\details{
  Usually, users will create distance specification using \code{\link{mdist}}, \code{\link{caliper}}, or
  \code{\link{exactMatch}}, but these functions are useful to generate sparse matching problems directly.
  If the data are already in a matrix form, use \code{as.InfinitySparseMatrix}. If you have the
  finite entries in a vector format, use \code{makeInfinitySparseMatrix}
}
\value{An object of class \code{InfinitySparseMatrix}, which will work as a distance argument
to \code{\link{fullmatch}} or \code{\link{pairmatch}}}
\author{Mark M. Fredrickson}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mdist}}, \code{\link{caliper}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{pairmatch}}
}
\examples{
  example.matrix <- matrix(c(1,2,Inf, 3,Inf,4, Inf,Inf,Inf), byrow = TRUE, nrow = 3,
                           dimnames = list(letters[1:3], LETTERS[24:26]))

  as.InfinitySparseMatrix(example.matrix)

  # create the same sparse matrix directly, function will create the appropriate dims
  # the data are in a different order, but the indices are correct
  (example.ism <- makeInfinitySparseMatrix(c(1,2,3,4), c(1,2,1,3), c(1,1,2,2),  LETTERS[24:26], letters[1:3]))
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
