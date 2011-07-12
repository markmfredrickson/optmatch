### maxControlsCap
require('optmatch')
data(nuclearplants)
mhd2a <- mdist(pr ~ date + cum.n, data = nuclearplants, 
              exclusions = exactMatch(pr ~ pt, data = nuclearplants))
mhd2a <- t(mhd2a)

mhd2a.caliper <- mhd2a + caliper(3, mhd2a)
stratumStructure(fullmatch(mhd2a.caliper)) # Works OK: 
maxControlsCap(mhd2a.caliper)              # no unmatchable Tx
stratumStructure(fullmatch(mhd2a.caliper, max=1))
stratumStructure(fullmatch(mhd2a.caliper), max=1/2)
stratumStructure(fullmatch(mhd2a.caliper)) # Problem in version <= .5-9:
(mx2 <- maxControlsCap(mhd2a.caliper))     # caused by unmatchable Tx
stratumStructure(fullmatch(mhd2a.caliper, max=mx2$strictest))
stratumStructure(fullmatch(mhd2a.caliper, max=1/2))
