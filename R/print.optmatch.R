print.optmatch <- function(x, quote=FALSE, ...)
{
   if (length(x) <= 0)
      cat("factor(0)\n")
   else 
   {
   rv <- as.character(x)
   names(rv) <- names(x)
   print(rv, quote=FALSE, ...)
   }
   invisible(x)
}
