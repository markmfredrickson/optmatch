require('optmatch')
data(plantdist)
# this will give a warning about not having data order
fullmatch(1 * (plantdist < 10)) # make plantdist < 10 numeric, not logical

data(nuclearplants)
mhd2 <- mdist(pr ~ date + cum.n, data = nuclearplants, 
              within = exactMatch(pr ~ pt, data = nuclearplants))
# the previous version of optmatch used fullmatch(mhd2 < 1)
# this is the equivalent using an ISM (logical operators treat them as numeric
# vectors)
mhd2[mhd2 < 1] <- 1
mhd2[mhd2 >= 1] <- 0
fullmatch(mhd2, data = nuclearplants)

