library(optmatch)
data(plantdist)

plantsfm <- fullmatch(plantdist) 
plantsfm1 <- fullmatch(plantdist,min.controls=2, max.controls=3)

stratumStructure(plantsfm)
stratumStructure(plantsfm1)

data(nuclearplants)
psd <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                       data=nuclearplants),
             standardization.scale = sd) # sd was standard in < 0.7-2
stratumStructure(fullmatch(psd))
stratumStructure(fm.psd.cal <- fullmatch(psd/(psd<.25)))

psd2 <- match_on(glm(pr~.-(pr+cost), family=binomial(),
                       data=nuclearplants),
              standardization.scale = sd,
              within = exactMatch(pr ~ pt, nuclearplants))
stratumStructure(fullmatch(psd2,min.controls=1,max.controls=1))
stratumStructure(fullmatch(psd2,min.controls=3,max.controls=3))


### Tests of min.controls, max.controls
stratumStructure(fm.psd.cal, min.controls=.5)
stratumStructure(fm.psd.cal, max.controls=3)
stratumStructure(fm.psd.cal, min.controls=.5, max.controls=3)

### tests of stratumStructure on non-optmatch objects
fac <- as.factor(fm.psd.cal)
tvec <- attr(fm.psd.cal, "contrast.group")
stratumStructure(fac, tvec)
stratumStructure(as.integer(fac),tvec)
