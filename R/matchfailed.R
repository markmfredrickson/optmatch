matchfailed <- function(matchobject) {
ans <- logical(length(matchobject))
if (!inherits(matchobject,'optmatch')) {
  warning('argument of matchfailed() not of class \'optmatch\'', call.=FALSE)
}
ans[grep("[.]NA$", as.character(matchobject))] <- TRUE
ans
}
