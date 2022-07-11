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
mhd <- match_on(pr ~ t1 + t2, data = nuclearplants) + caliper(match_on(ppty), width = 1)
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

## some random propensity scores
## set.seed(3943522034)
## n1 <- 30
## n0 <- 50
## z <- c(rep(1, n1), rep(0, n0))
## e <- rnorm(n1 + n0, sd = 0.5)
## p <- exp(z + e) / (1 + exp(z + e))
## pscores <- data.frame(z, round(p, 3))
pscores <- structure(list(z = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0), p = c(0.784, 0.777, 0.621, 0.744, 0.772,
                                                          0.82, 0.747, 0.72, 0.772, 0.834, 0.726, 0.738, 0.723, 0.773,
                                                          0.512, 0.835, 0.584, 0.664, 0.794, 0.695, 0.727, 0.742, 0.768,
                                                          0.724, 0.693, 0.856, 0.833, 0.684, 0.634, 0.549, 0.515, 0.453,
                                                          0.566, 0.336, 0.485, 0.361, 0.594, 0.451, 0.629, 0.64, 0.488,
                                                          0.54, 0.484, 0.516, 0.501, 0.644, 0.404, 0.388, 0.508, 0.578,
                                                          0.499, 0.445, 0.539, 0.512, 0.584, 0.671, 0.297, 0.45, 0.802,
                                                          0.36, 0.38, 0.515, 0.599, 0.406, 0.453, 0.471, 0.629, 0.597,
                                                          0.37, 0.524, 0.66, 0.307, 0.439, 0.614, 0.597, 0.487, 0.499,
                                                          0.525, 0.517, 0.433)), class = "data.frame", row.names = c(NA,
                                                                                                                     -80L))

## there is not much overlap of the propensity scores
boxplot(p ~ z, data = pscores)

## consequently the optimal match has many treated units sharing controls
m6 <- fullmatch(z ~ p, data = pscores)
stratumStructure(m6)

## we could force each set to include a maximum number of controls
## but another way to control the set sizes is to limit the maximum
## difference between the number of treated units sharing a control
## and the the number of sets containing more than one treated unit
m7 <- fullmatch(z ~ p, data = pscores, shared_treatment_excess = 9)
stratumStructure(m7)

## compare this to limiting sets to contain no more than 10 treated units
m8 <- fullmatch(z ~ p, data = pscores, min.controls = 1/10)
stratumStructure(m8)
