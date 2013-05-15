### scores is a wrapper for predict(), but if no newdata is specified, rather than using
### the original data, it looks for a data set in its parent environment, so should be
### called in a model.)

scores <- function(form, newdata=NULL,...)
{
  # If user didn't give newdata, extract from model call
  if(is.null(newdata))
    {
      # Get model call environment
      e <- parent.frame()
      # Unfortunately, it's basically called "attach(newdata)" - so lets re-merge it
      newdata <- data.frame(matrix(unlist(cbind(mget(ls(e),e))), nrow=length(get(ls(e)[1],e)), byrow=FALSE))
      # Note: Columns will be in alphabetical order now

      # Save the original names
      names(newdata) <- ls(e)
    }
  predict(form, newdata=newdata,...)
}
