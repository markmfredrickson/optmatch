### maxControlsCap
require('optmatch')
data(nuclearplants)
mhd2a <- mdist(pr ~ date + cum.n, data = nuclearplants, 
              excludes = exactMatch(pr ~ pt, data = nuclearplants))
mhd2a <- t(mhd2a)

stratumStructure(fullmatch(mhd2a/(mhd2a<3))) # Works OK: 
maxControlsCap(mhd2a/(mhd2a<3))              # no unmatchable Tx
stratumStructure(fullmatch(mhd2a/(mhd2a<3), max=1))
stratumStructure(fullmatch(mhd2a/(mhd2a<3), max=1/2))
stratumStructure(fullmatch(mhd2a/(mhd2a<2))) # Problem in version <= .5-9:
(mx2 <- maxControlsCap(mhd2a/(mhd2a<2)))     # caused by unmatchable Tx
stratumStructure(fullmatch(mhd2a/(mhd2a<2), max=mx2$strictest))
stratumStructure(fullmatch(mhd2a/(mhd2a<2), max=1/2))
