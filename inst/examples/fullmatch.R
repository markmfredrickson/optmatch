data(nuclearplants)
### Full matching on a Mahalanobis distance
mhd  <- match_on(pr ~ t1 + t2, data = nuclearplants)
( fm1 <- fullmatch(mhd, data = nuclearplants) )
summary(fm1)
### Full matching with restrictions
( fm2 <- fullmatch(mhd, min=.5, max=4, data = nuclearplants) )
summary(fm2)
### Full matching to half of available controls
( fm3 <- fullmatch(mhd, omit.fraction=.5, data = nuclearplants) )
summary(fm3)
### Full matching within a propensity score caliper.
ppty <- glm(pr~.-(pr+cost), family=binomial(), data = nuclearplants)
### Note that units without counterparts within the
### caliper are automatically dropped.
( fm4 <- fullmatch(mhd + caliper(match_on(ppty), 1), data = nuclearplants) )
summary(fm4)

### Propensity balance assessment. Requires RItools package.
if (require(RItools)) summary(fm4,ppty) 

### The order of the names in the match factor is the same
### as the nuclearplants data.frame since we used the data argument
### when calling fullmatch. The order would be unspecified otherwise.
cbind(nuclearplants, matches = fm4)
