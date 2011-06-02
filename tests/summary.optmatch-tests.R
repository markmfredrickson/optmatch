require('optmatch')
data(plantdist)
summary(fullmatch(plantdist<10))
summary(pairmatch(plantdist/(plantdist<1))) # Matching fails everywhere
data(nuclearplants)
psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
psd <- pscore.dist(psm)
psfm <- fullmatch(psd/(psd<.25))
summary(psfm)
pspm <- pairmatch(caliper(mdist(psm, structure.fmla=~pt), width=2)) # Fails in subclass '1'
summary(pspm)
psd$m[1,] <- psd$m[1,] + rep(100,22)
summary(pairmatch(psd, controls=2))
summary(psfm, propensity.model=psm)
require('RItools')
summary(psfm, propensity.model='foo')
summary(psfm, propensity.model=psm)
summary(psfm, psm)
psm2 <- glm(pr~ cut(date, c(67, 69.5, 72)) +
            t1 + t2 + cap + ne + ct + bw + cum.n + pt,
            family=binomial, data=nuclearplants)
psd2 <- pscore.dist(psm2)
summary(pairmatch(psd2), propensity.model=psm2)
summary(pspm, propensity.model=psm) # balance checking when matching has failed
                                    # in some subclasses
