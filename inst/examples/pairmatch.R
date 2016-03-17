data(nuclearplants)

### Pair matching on a Mahalanobis distance
( pm1 <- pairmatch(pr ~ t1 + t2, data = nuclearplants) )
summary(pm1)

### Pair matching within a propensity score caliper.
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
### For more complicated models, create a distance matrix and pass it to fullmatch.
mhd <- match_on(pr ~ t1 + t2, data = nuclearplants) + caliper(match_on(ppty), 2)
( pm2 <- pairmatch(mhd, data = nuclearplants) )
summary(pm2)

### Propensity balance assessment. Requires RItools package.
if(require(RItools)) summary(pm2, ppty)

### 1:2 matched triples
( tm <- pairmatch(pr ~ t1 + t2, controls = 2, data = nuclearplants) )
summary(tm)

### Creating a data frame with the matched sets attached.
### match_on(), caliper() and the like cooperate with pairmatch()
### to make sure observations are in the proper order:
all.equal(names(tm), row.names(nuclearplants))
### So our data frame including the matched sets is just
cbind(nuclearplants, matches=tm)

### In contrast, if your matching distance is an ordinary matrix
### (as earlier versions of optmatch required), you'll
### have to align it by observation name with your data set.
cbind(nuclearplants, matches = tm[row.names(nuclearplants)])


### Match in subgroups only. There are a few ways to specify this.
m1 <- pairmatch(pr ~ t1 + t2, data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m2 <- pairmatch(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
### Matching on propensity scores within matching in subgroups only:
m3 <- pairmatch(glm(pr ~ t1 + t2, data=nuclearplants, family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- pairmatch(glm(pr ~ t1 + t2 + pt, data=nuclearplants,
                    family=binomial),
                data=nuclearplants,
                within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- pairmatch(glm(pr ~ t1 + t2 + strata(pt), data=nuclearplants,
                    family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
