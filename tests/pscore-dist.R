require('optmatch')
data(nuclear, package="boot")
psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclear)
psd <- pscore.dist(psm)
fullmatch(psd/(psd<.25))
psm <- glm(pr~.-(pr+cost), family=binomial(), 
           data=nuclear, subset=(pt==0))      # v0.3-2 had trouble w/this
### fix was to change
###  makedist(structure.fmla,
###         data.frame(ZzZz, Ppty,glmobject$data),
### to
###  makedist(structure.fmla,
###         data.frame(ZzZz, Ppty,glmobject$data[names(Ppty),]),
psd0 <- pscore.dist(psm)
fullmatch(psd0/(psd0<.25))
psd1 <- pscore.dist(psm,standardization.scale=mad)
fullmatch(psd1/(psd1<.25))
