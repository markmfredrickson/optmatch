data(nuclearplants)
### Full matching on a Mahalanobis distance.
( fm1 <- fullmatch(pr ~ t1 + t2, data = nuclearplants) )
summary(fm1)

### Full matching with restrictions.
( fm2 <- fullmatch(pr ~ t1 + t2, min.controls = .5, max.controls = 4, data = nuclearplants) )
summary(fm2)

### Full matching to half of available controls.
( fm3 <- fullmatch(pr ~ t1 + t2, omit.fraction = .5, data = nuclearplants) )
summary(fm3)

### Full matching attempts recovery when the initial restrictions are infeasible.
### Limiting max.controls = 1 allows use of only 10 of 22 controls.
( fm4 <- fullmatch(pr ~ t1 + t2, max.controls = 1, data=nuclearplants) )
summary(fm4)
### To recover restrictions
optmatch_restrictions(fm4)

### Full matching within a propensity score caliper.
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
### Note that units without counterparts within the caliper are automatically dropped.
### For more complicated models, create a distance matrix and pass it to fullmatch.
mhd <- match_on(pr ~ t1 + t2, data = nuclearplants) + caliper(match_on(ppty))
( fm5 <- fullmatch(mhd, data = nuclearplants) )
summary(fm5)

### Propensity balance assessment. Requires RItools package.
if (require(RItools)) summary(fm5,ppty)

### The order of the names in the match factor is the same
### as the nuclearplants data.frame since we used the data argument
### when calling fullmatch. The order would be unspecified otherwise.
cbind(nuclearplants, matches = fm5)

### Match in subgroups only. There are a few ways to specify this.
m1 <- fullmatch(pr ~ t1 + t2, data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m2 <- fullmatch(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
### Matching on propensity scores within matching in subgroups only:
m3 <- fullmatch(glm(pr ~ t1 + t2, data=nuclearplants, family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- fullmatch(glm(pr ~ t1 + t2 + pt, data=nuclearplants,
                    family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- fullmatch(glm(pr ~ t1 + t2 + strata(pt), data=nuclearplants,
                    family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
