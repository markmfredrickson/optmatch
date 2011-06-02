require('optmatch')
data(plantdist)
pairmatch(plantdist)
pairmatch(plantdist, controls=2)
pairmatch(plantdist/(plantdist<1)) # Matching fails everywhere
pairmatch(plantdist/(plantdist<5), remove.unmatchables=TRUE) # Matching works after removing plant 'F'
data(nuclearplants)
psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
psd <- pscore.dist(psm)
pairmatch(psd, controls=2)
pairmatch(caliper(mdist(psm, structure.fmla=~pt), width=2)) # Fails in subclass '1'

