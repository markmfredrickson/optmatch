#' @export
#' @importFrom graphics boxplot
boxplot.glm <- function(x, data=NULL, xlab="Group", ylab=expression(paste(X, symbol("\242"), hat(beta))), main="Overlap on fitted scores",varwidth=TRUE, horizontal=FALSE, ...)
{
    # NB: if default xlab or ylab is altered,
    # update accordingly w/in function body below
    if (is.null(data))
      {
        dependent.variable <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
        linear.score <- x$linear.predictors
} else {
  linear.score <- predict(x, data, type = 'link', se.fit = FALSE)
  Data <- model.frame(terms(x), data)
  dependent.variable <- as.numeric(model.response(Data))
}
    if (horizontal) { #switch default axis labelings
        if (missing(xlab)) {
            xlab <- # default value of ylab
                expression(paste(X, symbol("\242"), hat(beta)))
        }
        if (missing(ylab)) {
            ylab <- # default value of xlab
                "Group"
        }
    }
boxplot(linear.score ~ dependent.variable, xlab=xlab, ylab=ylab,main=main, varwidth=varwidth,horizontal=horizontal,...)
  }

#' @export
boxplot.bigglm <- function(x, data,xlab="Group", ylab=expression(paste(X, symbol("\242"), hat(beta))), main="Overlap on fitted scores",varwidth=TRUE, horizontal=FALSE,...)
  {
  if (is.null(data)) {
    stop("data argument is required for making boxplots from bigglms")
  }

  if (!is.data.frame(data)) {
    stop("boxplot.bigglm doesn't understand data arguments that aren't data frames.")
  }

 linear.score <- predict(x, data, type = 'link', se.fit = FALSE)

  if (length(linear.score) != dim(data)[1]) {
    stop("predict.bigglm() returns a vector of the wrong length;
are there missing values in data?")
  }

  # this makes heavy use of the bigglm terms object, the original formula
  # if this implementation detail changes in later versions of bigglm,
  # look here for the problem.

  Data <-  model.frame(x$terms, data = data)
  dependent.variable <- as.numeric(model.response(Data))

  if (horizontal) { #switch default axis labelings
      if (missing(xlab)) {
          xlab <- # default value of ylab
              expression(paste(X, symbol("\242"), hat(beta)))
      }
      if (missing(ylab)) {
          ylab <- # default value of xlab
              "Group"
      }
  }
  boxplot(linear.score ~ dependent.variable, xlab=xlab, ylab=ylab,main=main, varwidth=varwidth, horizontal=horizontal,...)
  }
#' @export
#' @importFrom graphics boxplot
boxplot.svyglm <- function(x, xlab="Group", ylab=expression(paste(X, symbol("\242"), hat(beta))), main="Overlap on fitted scores",varwidth=TRUE, horizontal=FALSE, ...)
{
  if (requireNamespace("survey", quietly = TRUE)) {
    dependent.variable <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    dependent.variable  <- as.factor(dependent.variable)
        linear.score <- x$linear.predictors
    if (horizontal) { #switch default axis labelings
      if (missing(xlab)) {
          xlab <- # default value of ylab
              expression(paste(X, symbol("\242"), hat(beta)))
      }
      if (missing(ylab)) {
          ylab <- # default value of xlab
              "Group"
      }
    }
    survey::svyboxplot(linear.score ~ dependent.variable, design=x$survey.design, main=main, varwidth=varwidth, horizontal=horizontal, ...)
  } else {
    stop("'survey' package must be installed")
  }
  NULL
}
