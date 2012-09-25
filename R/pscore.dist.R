pscore.dist <- function(glmobject, structure.fmla=NULL,standardization.scale=sd)
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
