### maxControlsCap
require('optmatch')
data(nuclearplants)
mhd2a <- match_on(pr ~ date + cum.n, data = nuclearplants, 
              exclusions = exactMatch(pr ~ pt, data = nuclearplants))
mhd2a <- t(mhd2a)

mhd2a.caliper <- mhd2a + caliper(mhd2a, 3)
stratumStructure(fullmatch(mhd2a.caliper)) # Works OK: 
maxControlsCap(mhd2a.caliper)              # no unmatchable Tx
stratumStructure(fullmatch(mhd2a.caliper, max=1))
stratumStructure(fullmatch(mhd2a.caliper, max=1/2))
stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2))) # Problem in version <= .5-9:
(mx2 <- maxControlsCap(mhd2a + caliper(mhd2a, 2)))     # caused by unmatchable Tx
stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=mx2$strictest))
stratumStructure(fullmatch(mhd2a + caliper(mhd2a, 2), max=1/2))
