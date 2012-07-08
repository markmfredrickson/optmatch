mahal.dist <- function(distance.fmla, data, structure.fmla=NULL, inverse.cov=NULL)
  {
    
if (is.null(structure.fmla))
  {
  if (!attr(terms(distance.fmla,data=data), "response")>0) 
    stop("either distance.fmla or structure.fmla must specify a treatment group variable")
  structure.fmla <- update.formula(distance.fmla, .~1,data=data)
  structure.fmla <- terms.formula(structure.fmla, data=data)
} else
{
    if (!attr(terms(structure.fmla,data=data), "response")>0 &
        !attr(terms(distance.fmla,data=data), "response")>0)
      stop("either distance.fmla or structure.fmla must specify a treatment group variable")

    if (!attr(terms(structure.fmla,data=data), "response")>0)
      {
      lhs <- as.character(distance.fmla[[2]])
      structure.fmla <- update.formula(structure.fmla,
                                       paste(lhs, '~.', sep=''))
      structure.fmla <- terms.formula(structure.fmla, data=data)
      }
  }
distance.fmla <- update.formula(distance.fmla, ~-1+.,data=data)
distance.fmla <- terms(distance.fmla, data=data)
ds.vars <-all.vars(distance.fmla)
if (length(ds.vars)<2) stop("No variables on RHS of distance.fmla")

inp <- parse(text=paste("list(", paste(ds.vars, collapse = ","), ")")) # from def of get_all_vars()

dfr <- model.matrix(distance.fmla, #model.frame(distance.fmla,data))
                    structure(eval(inp, data, parent.frame()),
                              names=as.character(ds.vars)))
sf.vars <- all.vars(structure.fmla)
sf.vars <- sf.vars[!(sf.vars%in% names(dfr))]
sf.vars <- sf.vars[sf.vars %in% names(data)]

if (is.null(inverse.cov))
  {
    zpos <- attr(terms(structure.fmla,data=data), "response")
    vars <- eval(attr(terms(structure.fmla,data=data), "variables"), data, 
                 parent.frame())
    zz <- vars[[zpos]]
    if (!(is.numeric(zz) || is.logical(zz)))
      stop("Treatment variable should be logical or numeric")
    zz <- zz > 0
    
  cv <- cov(dfr[as.logical(zz), ,drop=FALSE])*(sum(zz)-1)/(length(zz)-2)
  cv <- cv + cov(dfr[!zz,,drop=FALSE])*(sum(!zz)-1)/(length(zz)-2)
  icv <- try( solve(cv), silent=TRUE)
  if (inherits(icv,"try-error"))
    {
       dnx <- dimnames(cv)
       s <- svd(cv)
       nz <- (s$d > sqrt(.Machine$double.eps) * s$d[1])
       if (!any(nz))
         stop("covariance has rank zero")

       icv <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
       dimnames(icv) <- dnx[2:1]
    }
## stopifnot(all.equal(dimnames(icv)[[1]],dimnames(dfr)[[2]]))
  } else
{
  if (!is.matrix(inverse.cov) || dim(inverse.cov)[1]!=dim(inverse.cov)[2] ||
      (!is.null(dimnames(inverse.cov)) &&
       !isTRUE(all.equal(dimnames(inverse.cov)[[1]],dimnames(inverse.cov)[[2]]))
       )
      ) stop("inverse.cov must be a square symmetric matrix.")

  if (dim(inverse.cov)[1]!=dim(dfr)[2]) stop("dimension of inverse.cov must match number of data terms.")
  if (!is.null(dimnames(inverse.cov)) &&
      !isTRUE(all.equal(dimnames(inverse.cov)[[1]],dimnames(dfr)[[2]]))
      )
    {
      icvnm <- gsub("TRUE$", "", dimnames(inverse.cov)[[1]])
      dfrnm <- gsub("TRUE$", "", dimnames(dfr)[[2]])
      if (!isTRUE(all.equal(icvnm, dfrnm)))    stop("dimnames of inverse.cov don't match names of data terms") 
    }
  if (is.null(dimnames(inverse.cov)) & dim(inverse.cov)[1] > 1) warning("inverse.cov lacks dimnames so I can't confirm it's aligned with data terms.")
  
icv <- inverse.cov
}
attr(structure.fmla, 'generation.increment') <- 1

ln.dfr <- dim(dfr)[2]
dfr <- data.frame(dfr, data[sf.vars], row.names=row.names(data))
dimnames(icv) <- list(names(dfr)[1:ln.dfr], names(dfr)[1:ln.dfr])

makedist(structure.fmla, dfr,
         fn=optmatch.mahalanobis, inverse.cov=icv)
  }



optmatch.mahalanobis <- function(trtvar, dat, inverse.cov)
  {  
  myMH <- function(Tnms, Cnms, inv.cov, data) {
   stopifnot(!is.null(dimnames(inv.cov)[[1]]), 
             all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
             all(dimnames(inv.cov)[[1]] %in% names(data)))
   covars <- dimnames(inv.cov)[[1]]
   xdiffs <- as.matrix(data[Tnms,covars])
   xdiffs <- xdiffs - as.matrix(data[Cnms,covars])
   rowSums((xdiffs %*% inv.cov) * xdiffs)
 }

te <- try(ans <- outer(names(trtvar)[trtvar], names(trtvar)[!trtvar],
               FUN=myMH, inv.cov=inverse.cov, data=dat), silent=TRUE)

if (inherits(te, 'try-error'))
  {
    if (substr(unclass(te),1,29)!="Error: cannot allocate vector")
      stop(unclass(te))
ans <- matrix(0,sum(trtvar), sum(!trtvar))

  cblocks <- 1
while (inherits(te, 'try-error') &&
       substr(unclass(te),1,29)=="Error: cannot allocate vector" &&
       (sum(!trtvar)/cblocks)>1)
  {
cblocks <- cblocks*2
bsz <- ceiling(sum(!trtvar)/cblocks)
trtvar.ctlnms <-
  split(names(trtvar)[!trtvar],
        rep(1:cblocks,rep(bsz, cblocks))[1:sum(!trtvar)]
        )
te <- try(lapply(trtvar.ctlnms,
                 function(ynms) outer(names(trtvar)[trtvar],ynms,
                                      FUN=myMH, inv.cov=inverse.cov,
                                      data=dat) ), silent=TRUE )
###for (ii in seq(sum(!trtvar), 1, by=bsz) )
###  {
###    ans[,max((ii-bsz+1),1):ii] <-
###      outer(names(trtvar)[trtvar],
###           trtvar.ctlnms[max((ii-bsz+1),1):ii],
###            FUN=myMH, inv.cov=inverse.cov, data=dat)
###  }, silent=TRUE)
}
if (inherits(te, 'try-error') )
      {stop(unclass(te)) } else ans <- unlist(te)
          
  }
dim(ans) <- c(sum(trtvar), sum(!trtvar))

  dimnames(ans) <- list(names(trtvar)[trtvar], names(trtvar)[!trtvar])
  ans
  }
