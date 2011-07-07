require('optmatch')
data(plantdist)
fullmatch(plantdist<10)
data(nuclearplants)
mhd2 <- mdist(pr ~ date + cum.n, data = nuclearplants, 
              excludes = exactMatch(pr ~ pt, data = nuclearplants))
fullmatch(mhd2 < 1)
