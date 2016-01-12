setOldClass(c("optmatch.dlist", "list"))

##' Deprecated in favor of \code{\link{match_on}}
##'
##' The \code{mdist} method provides three ways to construct a
##' matching distance (i.e., a distance matrix or suitably organized
##' list of such matrices): guided by a function, by a fitted model,
##' or by a formula.  The class of the first argument given to
##' \code{mdist} determines which of these methods is invoked.
##'
##' The \code{mdist.function} method takes a function of two
##' arguments. When called, this function will receive the treatment
##' observations as the first argument and the control observations as
##' the second argument. As an example, the following computes the raw
##' differences between values of \code{t1} for treatment units (here,
##' nuclear plants with \code{pr==1}) and controls (here, plants with
##' \code{pr==0}), returning the result as a distance matrix:
##'
##' \code{sdiffs <- function(treatments, controls) {
##'      abs(outer(treatments$t1, controls$t1, `-`))
##'    }
##'  }
##'
##'  The \code{mdist.function} method does similar things as the
##'  earlier optmatch function \code{makedist}, although the interface
##'  is a bit different.
##'
##'  The \code{mdist.formula} method computes the squared Mahalanobis
##'  distance between observations, with the right-hand side of the
##'  formula determining which variables contribute to the Mahalanobis
##'  distance. If matching is to be done within strata, the
##'  stratification can be communicated using either the
##'  \code{structure.fmla} argument (e.g. \code{~ grp}) or as part of
##'  the main formula (e.g. \code{z ~ x1 + x2 | grp}).
##'
##'  An \code{mdist.glm} method takes an argument of class \code{glm}
##'  as first argument.  It assumes that this object is a fitted
##'  propensity model, extracting distances on the linear propensity
##'  score (logits of the estimated conditional probabilities) and, by
##'  default, rescaling the distances by the reciprocal of the pooled
##'  s.d. of treatment- and control-group propensity scores.  (The
##'  scaling uses \code{mad}, for resistance to outliers, by default;
##'  this can be changed to the actual s.d., or rescaling can be
##'  skipped entirely, by setting argument
##'  \code{standardization.scale} to \code{sd} or \code{NULL},
##'  respectively.)  A \code{mdist.bigglm} method works analogously
##'  with \code{bigglm} objects, created by the \code{bigglm} function
##'  from package \sQuote{biglm}, which can handle bigger data sets
##'  than the ordinary glm function can.  In contrast with
##'  \code{mdist.glm} it requires additional \code{data} and
##'  \code{structure.fmla} arguments.  (If you have enough data to
##'  have to use \code{bigglm}, then you'll probably have to subgroup
##'  before matching to avoid memory problems. So you'll have to use
##'  the \code{structure.fmla} argument anyway.)
##'
##' @title (Deprecated, in favor of \code{\link{match_on}}) Create
##'   matching distances
##' @param x The object to use as the basis for forming the mdist.
##'   Methods exist for formulas, functions, and generalized linear
##'   models.
##' @param structure.fmla A formula denoting the treatment variable on
##'   the left hand side and an optional grouping expression on the
##'   right hand side. For example, \code{z ~ 1} indicates no
##'   grouping. \code{z ~ s} subsets the data only computing distances
##'   within the subsets formed by \code{s}. See method notes, below,
##'   for additional formula options.
##' @param data Data where the variables references in `x` live.
##' @param subset If non-NULL, the subset of `data` to be used.
##' @param standardization.scale A function to scale the distances; by
##'   default uses `mad`.
##' @param trtgrp Dummy variable for treatment group membership.
##' @param ... Additional method arguments. Most methods require a
##'   'data' argument.
##' @return Object of class \code{optmatch.dlist}, which is suitable
##'   to be given as \code{distance} argument to
##'   \code{\link{fullmatch}} or \code{\link{pairmatch}}.
##' @author Mark M. Fredrickson
##' @references P.~R. Rosenbaum and D.~B. Rubin (1985),
##'   \sQuote{Constructing a control group using multivariate matched
##'   sampling methods that incorporate the propensity score},
##'   \emph{The American Statistician}, \bold{39} 33--38.
##' @seealso \code{\link{fullmatch}}, \code{\link{pairmatch}},
##'   \code{\link{match_on}}
##' @keywords nonparametric
##' @export
##' @rdname mdist
mdist <- function(x, structure.fmla = NULL, ...) {
  cl <- match.call()
  UseMethod("mdist", x)
}
getCall.optmatch.dlist <- function(x, ...) attr(x, "call")


##' @export
##' @rdname mdist
mdist.optmatch.dlist <- function(x, structure.fmla = NULL, ...) {
  return(x)
} # just return the argument

# mdist method: function
# for the function method, both data and structure fmla are required
# api change from makedist: I think it would make more sense to have
# the function take two data.frames: the treatments in this stratum
# and the controls in the stratum. It could then return the matrix of
# mdists, which the rest of the function would markup with rownames
# etc.
##' @export
##' @rdname mdist
mdist.function <- function(x, structure.fmla = NULL, data = NULL, ...) {

  if (is.null(data) || is.null(structure.fmla)) {
    stop("Both data and the structure formula are required for
    computing mdists from functions.")
  }
  if (!exists("cl")) cl <- match.call()
  theFun <- match.fun(x)
  parsedFmla <- parseFmla(structure.fmla)

  if(is.null(parsedFmla[[1]])) {
    stop("Treatment variable required")
  }

  if((identical(parsedFmla[[2]], 1) && is.null(parsedFmla[3])) ||
     (length(parsedFmla[[2]]) > 1)) {
    stop("Please specify the grouping as either: z ~ grp or z ~ 1 | grp")
  }

  treatmentvar <- parsedFmla[[1]]

  if(is.null(parsedFmla[3][[1]])) { # I swear subscripting is incomprehensible!
    strata <- parsedFmla[[2]]
  } else {
    strata <- parsedFmla[[3]]
  }


  # split up the dataset by parts
  # call optmatch.dlist maker function on parts, names

    # create a function to produce one distance matrix
  doit <- function(data,...) {
     # indicies are created per chunk to split out treatment and controls
     indices <- data[as.character(treatmentvar)] == 1
     treatments <- data[indices,]
     controls <- data[!indices,]
     distances <- theFun(treatments, controls, ...)

     colnames(distances) <- rownames(controls)
     rownames(distances) <- rownames(treatments)

     return(distances)
  }

  if (!(identical(strata, 1))) {
    if(is.factor(eval(strata, data))){
      ss <- eval(strata, data) ##to preserve existing labels/levels
    } else {
      ss <- factor(eval(strata, data), labels = 'm')
    }
    ans <- lapply(split(data, ss), doit,...)
  } else {
    ans <- list(m = doit(data,...))
  }

  attr(ans, 'row.names') <- row.names(data)
  attr(ans, "call") <- cl

  class(ans) <- c('optmatch.dlist', 'list')
  return(ans)
}

##' @export
##' @rdname mdist
mdist.formula <- function(x, structure.fmla = NULL, data = NULL, subset=NULL,...) {
  mf <- match.call(expand.dots=FALSE)
  if (!exists("cl")) cl <- match.call()
  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=="x"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")

  if (length(x)==2)
      if (!is.null(structure.fmla) && length(structure.fmla)==3)
          x <- update.formula(structure.fmla, x) else
  stop("If you give mdist() a formula, it needs to have a left hand side\nin order for mdist.formula() to figure out who is to be matched to whom.",
       call.=FALSE)

  if (isThereAPipe(x)) {
      if (!is.null(structure.fmla))
          warning("I see a pipe, a '|', in the formula you gave as primary argument to mdist().\n So I'm using what's to the right of it in lieu of the RHS of the\n'structure.fmla' argument that you've also given.", call.=FALSE)
      parsed <- parseFmla(x)
    # this block occurs if the grouping factor is present
    # e.g. z ~ x1 + x2 | grp
    if (!is.null(unlist(parsed[3]))) {
      xenv <- environment(x)
      x <- as.formula(paste(as.character(parsed[1:2]), collapse = "~"))
      environment(x) <- xenv
      structure.fmla <- as.formula(paste("~", parsed[[3]]))
    }
  }
  mf$formula <-  makeJoinedFmla(x, structure.fmla)
  mf <- eval(mf, parent.frame())

###  return(mf)
  ans <- old.mahal.dist(x, data = mf, structure.fmla = structure.fmla, ...)
  attr(ans, "call") <- cl
  ans
}

isThereAPipe <- function(fmla)
{
    inherits(fmla, "formula") && length(fmla) == 3 && length(fmla[[3]])==3 && fmla[[3]][[1]] == as.name("|")
}
# One big formula from a formula plus astructure formula -- for use w/ model.frame
makeJoinedFmla <- function(fmla, structure.fmla)
{
    if (is.null(structure.fmla)) return(fmla)
    stopifnot(inherits(fmla, "formula"), inherits(structure.fmla, "formula"),
              length(fmla)==3)

l <- length(structure.fmla)
structure.fmla[[l]] <- as.call(c(as.name("+"), as.name("."), structure.fmla[[l]]))

update.formula(fmla, structure.fmla)
}

##' @export
##' @rdname mdist
mdist.glm <- function(x, structure.fmla = NULL, standardization.scale=mad, ...)
{
  if (!exists("cl")) cl <- match.call()
  ans <- old.pscore.dist(x,  structure.fmla = structure.fmla, standardization.scale=standardization.scale, ...)
  attr(ans, "call") <- cl
  ans
}

old.pscore.dist <- function(glmobject, structure.fmla=NULL,standardization.scale=sd)
  {
stopifnot(all(c('y', 'linear.predictors','data')%in%names(glmobject)))

if (is.null(structure.fmla))
  {
    structure.fmla <- as.formula("ZzZz~1")
} else
{
structure.fmla <- update.formula(structure.fmla, ZzZz~.)
if (!all(all.vars(structure.fmla)%in%c('ZzZz',names(model.frame(glmobject)))))
  warning('stratifying variables (in structure.fmla) not in propensity specification')
}
ZzZz <- glmobject$y>0
pooled.sd <- if (is.null(standardization.scale)){
  1 } else szn.scale(glmobject$linear.predictors,ZzZz,standardization.scale)
PpTy <- glmobject$linear.predictors/pooled.sd

attr(structure.fmla, 'generation.increment') <- 1

makedistOptmatchDlist(structure.fmla,
         data.frame(ZzZz, PpTy,model.frame(glmobject)),
         fn=function(trtvar,data)
         {
           sclr <- data[names(trtvar), 'PpTy']
           names(sclr) <- names(trtvar)
           abs(outer(sclr[trtvar], sclr[!trtvar], '-'))
         }
           )

}

szn.scale <- function(x,Tx,standardizer=mad,...) {
sqrt( ((sum(!Tx)-1)*standardizer(x[!Tx])^2 +
       (sum(!!Tx)-1)*standardizer(x[!!Tx])^2)/
     (length(x)-2)
     )
}

# parsing formulas for creating mdists
parseFmla <- function(fmla) {

  treatment <- fmla[[2]]
  rhs <- fmla[[3]]
  if (length(rhs) == 3 && rhs[[1]] == as.name("|")) {
    covar <- rhs[[2]]
    group <- rhs[[3]]

  } else {
    covar <- rhs
    group <- NULL
  }

  return(c(treatment, covar, group))

}

##' @export
##' @rdname mdist
mdist.bigglm <- function(x, structure.fmla = NULL, data = NULL, standardization.scale=mad, ...)
{
  if (is.null(data))
    stop("data argument is required for computing mdists from bigglms")

  if (!is.data.frame(data))
    stop("mdist doesn't understand data arguments that aren't data frames.")

  if (is.null(structure.fmla) | !inherits(structure.fmla, 'formula'))
    stop("structure.fmla argument required with bigglms.
(Use form 'structure.fmla=<treatment.variable> ~ 1'
 for no stratification before matching)")
  if (!exists("cl")) cl <- match.call()

theps <- predict(x, data, type='link', se.fit=FALSE)
if (length(theps)!=dim(data)[1])
stop("predict.bigglm() returns a vector of the wrong length;
are there missing values in data?")


Data <-  model.frame(structure.fmla, data=data)
treatmentvar <- as.character(structure.fmla[[2]])
pooled.sd <- if (is.null(standardization.scale)) 1 else
szn.scale(theps, Data[[treatmentvar]], standardizer=standardization.scale,...)

Data$tHePs <- theps/pooled.sd

psdiffs <- function(treatments, controls) {
abs(outer(as.vector(treatments$tHePs),
as.vector(controls$tHePs), `-`))
}

ans <- mdist(psdiffs, structure.fmla=structure.fmla,
      data=Data)
  attr(ans, "call") <- cl
  ans
}


### mdist method: numeric.
### (mdist can't work with numeric vectors at present,
### but it can return an informative error message).
##' @export
##' @rdname mdist
mdist.numeric <- function(x, structure.fmla = NULL, trtgrp=NULL, ...) {
  stop("No method exists for 'numeric' objects. See 'mdist.formula' for an alternative.")
}

makedistOptmatchDlist <- function(structure.fmla, data,
                     fn=function(trtvar, dat, ...){
                       matrix(0, sum(trtvar), sum(!trtvar),
                              dimnames=list(names(trtvar)[trtvar],
                                names(trtvar)[!trtvar]))},
                     ...)
  {
  if (!attr(terms(structure.fmla), "response")>0)
    stop("structure.fmla must specify a treatment group variable")
  fn <- match.fun(fn)

### WHEN THIS FUNCTION IS WRAPPED TO, THIS IS HOW INFO ABOUT WHICH
### GENERATION PARENT FRAME structure.fmla IS TO BE EVALUATED IN IS PASSED
  pframe.generation <- 1
  if (!is.null(attr(structure.fmla, "generation.increment")))
    pframe.generation <- pframe.generation +
      attr(structure.fmla, "generation.increment")

  zpos <- attr(terms(structure.fmla), "response")
  vars <- eval(attr(terms(structure.fmla), "variables"), data,
             parent.frame(n=pframe.generation))
  zzz <- vars[[zpos]]
  if (!is.numeric(zzz) & !is.logical(zzz))
    stop("treatment variable (LHS of structure.fmla) must be numeric or logical")
  if (any(is.na(zzz)))
      stop("NAs not allowed in treatment variable (LHS of structure.fmla)")
  if (all(zzz>0))
    stop("there are no controls (LHS of structure.fmla >0)")
  if (all(zzz<=0))
    stop("there are no treatment group members (LHS of structure.fmla <=0)")

  zzz <- (zzz>0)
  vars <- vars[-zpos]
  names(zzz) <- row.names(data)
if (length(vars)>0)
  {
  ss <- interaction(vars, drop=TRUE)
} else ss <- factor(zzz>=0, labels="m")
  ans <- tapply(zzz, ss, FUN=fn,
                dat=data, ..., simplify=FALSE)
  FUNchk <- unlist(lapply(ans,
                          function(x){!is.matrix(x) & !is.vector(x)}))
  if (any(FUNchk)) { stop("fn should always return matrices")}

  mdms <- split(zzz,ss)
  NMFLG <- FALSE

  for (ii in (1:length(ans)))
  {
    dn1 <- names(mdms[[ii]])[mdms[[ii]]]
    dn2 <- names(mdms[[ii]])[!mdms[[ii]]]
    if (is.null(dim(ans[[ii]])))
      {
      if (length(dn1)>1 & length(dn2)>1)
        { stop("fn should always return matrices")}
      if (length(ans[[ii]])!=max(length(dn1), length(dn2)))
        { stop(paste("unuseable fn value for stratum", names(ans)[ii]))}

      if (is.null(names(ans[[ii]])))
          {
           ans[[ii]] <-  matrix(ans[[ii]], length(dn1), length(dn2),
                         dimnames=list(dn1,dn2))
         } else
        { if (length(dn1)>1)
            {
              ans[[ii]] <- matrix(ans[[ii]],length(dn1), 1,
                                 dimnames=list(names(ans[[ii]]), dn2))
            } else {
              ans[[ii]] <- matrix(ans[[ii]], 1, length(dn2),
                                 dimnames=list(dn1, names(ans[[ii]])))
            }
        }
      } else {
      if (!all(dim(ans[[ii]])==c(length(dn1), length(dn2))))
        { stop(paste("fn value has incorrect dimension at stratum",
                     names(ans)[ii])) }
      if (is.null(dimnames(ans[[ii]])))
        {
          dimnames(ans[[ii]]) <- list(dn1, dn2)
          NMFLG <- TRUE
        } else {
          if (!all(dn1%in%dimnames(ans[[ii]])[[1]]) |
              !all(dn2%in%dimnames(ans[[ii]])[[2]]) )
            { stop(paste(
                    "dimnames of fn value don't match unit names in stratum",
                         names(ans)[ii])) }
        }

      }
  }

  if (NMFLG){
warning("fn value not given dimnames; assuming they are list(names(trtvar)[trtvar], names(trtvar)[!trtvar])")}


  attr(ans, 'row.names') <- names(zzz)
  class(ans) <- c('optmatch.dlist', 'list')
  ans
  }

old.mahal.dist <- function(distance.fmla, data, structure.fmla=NULL, inverse.cov=NULL)
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

makedistOptmatchDlist(structure.fmla, dfr,
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
