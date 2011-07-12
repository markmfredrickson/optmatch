require('optmatch')
data(plantdist)
fullmatch(plantdist<10)
data(nuclearplants)
mhd2 <- mdist(pr ~ date + cum.n, data = nuclearplants, 
              exclusions = exactMatch(pr ~ pt, data = nuclearplants))
fullmatch(caliper(1, mhd2))
