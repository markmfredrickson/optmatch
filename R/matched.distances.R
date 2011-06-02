matched.distances <- function(matchobj, distance, preserve.unit.names=FALSE)
  {
stopifnot(inherits(matchobj,"optmatch"))

finddist.mat <- function(dmat, omobj)
  {
tapply(names(omobj),omobj, FUN=function(x,DMAT){
  DMAT[match(x,dimnames(DMAT)[[1]], nomatch=0),
                 match(x,dimnames(DMAT)[[2]], nomatch=0),
       drop=!preserve.unit.names]
}, dmat)
  }

if (!inherits(distance,"optmatch.dlist"))
  {
return(finddist.mat(distance, matchobj))
  } else {
res <- lapply(distance,finddist.mat, matchobj)
res <- lapply(res, function(x){x[unlist(lapply(x,length))>0]})
names(res) <- NULL
nms <- unlist(lapply(res, names))
if (any(duplicated(nms)))
warning("something is wrong -- matched set referenced in separate distance matrices")
res <- unlist(res, recursive=FALSE)
return(res[levels(matchobj)])
  }

  }
