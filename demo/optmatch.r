library(optmatch)

pause <- function(str) {
  invisible(readline(paste("\nNext:", str, "<press ret>"))) 
}

data(nuclearplants)

pause("A propensity score distance")
(dist.glm <- mdist(glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)))

pause("A Mahalanobis distance")
(dist.mahal <- mdist(pr ~ t1 + t2, data = nuclearplants))

pause("Absolute difference on a scalar-distance")

sdiffs <- function(treatments, controls) {
  abs(outer(treatments$t1, controls$t1, `-`))
}

(dist.abs <- mdist(sdiffs, structure.fmla = pr ~ 1, data = nuclearplants))

pause("Pair match")
plantspm <- pairmatch(dist.glm) 
stratumStructure(plantspm)

(plantspm.d <- matched.distances(plantspm,dist.glm))
mean(unlist(plantspm.d))

pause("A 1:2 match")
plantstm <- pairmatch(dist.glm,controls=2) 
stratumStructure(plantstm)
mean(unlist(matched.distances(plantstm,dist.glm)))
unlist(lapply(matched.distances(plantstm,dist.glm), max))

pause("A full match")
plantsfm <- fullmatch(dist.mahal) 
stratumStructure(plantsfm)
mean(unlist(matched.distances(plantsfm,dist.glm)))

pause("Full matching with restrictions")
plantsfm1 <- fullmatch(dist.mahal, 
min.controls=2, max.controls=3)   
stratumStructure(plantsfm1)
mean(unlist(matched.distances(plantsfm1,dist.mahal)))

all(plantsfm1[sort(names(plantsfm1))]
  == fullmatch(
       t(dist.mahal$m),
       min.controls=(1/3),
       max.controls=(1/2)) )


pause("Mahalanobis matching with propensity calipers")

plantsfm2 <- fullmatch(caliper(dist.glm, width = .25) + dist.mahal)

matched.distances(plantsfm2,dist.glm)
matched.distances(plantsfm2,dist.mahal)

pause("Combining matches with original data")

### The usual scenario is that rownames will be in the same
### order as the original data, though it safe to check:
all.equal(names(plantsfm2), row.names(nuclearplants))
cbind(nuclearplants, matches=plantsfm2)

### If the match is not in original order:
cbind(nuclearplants, matches = plantsfm2[row.names(nuclearplants)])

