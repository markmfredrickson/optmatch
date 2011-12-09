\name{mdist-methods}
\docType{methods}
\alias{mdist-methods}
\alias{mdist,bigglm-method}
\alias{mdist,formula-method}
\alias{mdist,function-method}
\alias{mdist,glm-method}
\alias{mdist,InfinitySparseMatrix-method}
\alias{mdist,matrix-method}
\alias{mdist,numeric-method}
\title{Methods for Function \code{mdist} in Package \pkg{optmatch}}
\description{Create
  matrices of distances between observations to be used in the match process.
  Using these matrices, 
  \code{pairmatch()} or \code{fullmatch()} can determine optimal matches.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "formula")}}{
The formula method produces, by default, a Mahalanobis distance specification based
on the formula \code{treatment ~ var1 + var2 + ... }, where \code{treatment} is the
treatment indicator. A Mahalanobis distance scales the squared Euclidean distance by
the inverse of the covariance matrix. Other scale matrices can be supplied, for example,
the identity matrix will result in squared Euclidean distance. If you wish to use an
alternative function to compute the covariance, pass it in the \code{COV} argument.
}

\item{\code{signature(x = "function")}}{
  Additionally requires the \code{z} (a treatment indicator) and \code{data} (a vector, matrix, or data.frame).
  The function takes two arguments, a vector representing the treated value for all treated-control distances and
  a vector of the same length representing the control value. It should return a vector of the same length representing
  the distances between each treated-control pair.
}

\item{\code{signature(x = "glm")}}{
  Produces a propensity score distance matrix. Standardizes the sum of the linear predictors
  by the supplied \code{standardization.scale} argument. By default uses median absolute deviation.
}

\item{\code{signature(x = "bigglm")}}{
  Same as \code{glm} method, but for \code{bigglm} objects.
}

\item{\code{signature(x = "InfinitySparseMatrix")}}{
  Returns the argument directly.
}

\item{\code{signature(x = "matrix")}}{
  Returns the argument directly.
}

\item{\code{signature(x = "numeric")}}{
  Only included to return a helpful error message and suggestions for producing
  a distance specification from a single numeric vector.
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
