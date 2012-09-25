abs.optmatch.dlist <- function (x) 
{
    
     rn1 <- attr(x, "row.names")
     x <- x[unlist(as.logical(lapply(x, length)))]

    value <- lapply(x, abs)

    class(value) <- c('optmatch.dlist', 'list')
        attr(value, "row.names") <- rn1
    value
  }
