require('optmatch')

data(plantdist)
# warnings about lacking data argument are expected
pairmatch(plantdist)
pairmatch(plantdist, controls=2)
## Not run: pairmatch(plantdist + caliper(plantdist, 1)) # Matching fails everywhere
pairmatch(plantdist + caliper(plantdist, 5, compare = `<`), remove.unmatchables=TRUE) # Matching works after removing plant 'F'

data(nuclearplants)
# in both of mdist calls below use sd to maintain backwards compatibility with
# pscore.dist, which used sd by default. mdist has used mad as the std. scale
# since it was added to the package, so the use of mdist should be consistent
# for users going forward.
psm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
psd <- mdist(psm, standardization.scale = sd)
pairmatch(psd, controls=2, data = nuclearplants)

# Also not run: again an error would be thrown (which R CMD CHECK does not like)
# pairmatch(caliper(mdist(psm, standardization.scale = sd,
#   within = exactMatch(pr ~ pt, data =
#   nuclearplants)), width=2)) # Fails in subclass '1'

