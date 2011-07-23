matched <- function(matchobject) {
ans <- logical(length(matchobject))
if (!inherits(matchobject,'optmatch')) {
  warning('argument of matched() not of class \'optmatch\'', call.=FALSE)
        }
ans[grep("[.][123456789][0123456789]*$", as.character(matchobject))] <- TRUE
ans
}
