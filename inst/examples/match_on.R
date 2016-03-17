data(nuclearplants)
match_on.examples <- list()
### Propensity score distances.
### Recommended approach:
(aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants))
match_on.examples$ps1 <- match_on(aGlm)
### A second approach: first extract propensity scores, then separately
### create a distance from them.  (Useful when importing propensity
### scores from an external program.)
plantsPS <- predict(aGlm)
match_on.examples$ps2 <- match_on(pr~plantsPS, data=nuclearplants)
### Full matching on the propensity score.
fm1 <- fullmatch(match_on.examples$ps1, data = nuclearplants)
fm2 <- fullmatch(match_on.examples$ps2, data = nuclearplants)
### Because match_on.glm uses robust estimates of spread,
### the results differ in detail -- but they are close enough
### to yield similar optimal matches.
all(fm1 == fm2) # The same

### Mahalanobis distance:
match_on.examples$mh1 <- match_on(pr ~ t1 + t2, data = nuclearplants)

### Absolute differences on a scalar:
tmp <- nuclearplants$t1
names(tmp) <- rownames(nuclearplants)

(absdist <- match_on(tmp, z = nuclearplants$pr,
                  within = exactMatch(pr ~ pt, nuclearplants)))

### Pair matching on the variable `t1`:
pairmatch(absdist, data = nuclearplants)


### Propensity score matching within subgroups:
match_on.examples$ps3 <- match_on(aGlm, exactMatch(pr ~ pt, nuclearplants))
fullmatch(match_on.examples$ps3, data = nuclearplants)

### Propensity score matching with a propensity score caliper:
match_on.examples$pscal <- match_on.examples$ps1 + caliper(match_on.examples$ps1, 1)
fullmatch(match_on.examples$pscal, data = nuclearplants) # Note that the caliper excludes some units

### A Mahalanobis distance for matching within subgroups:
match_on.examples$mh2 <- match_on(pr ~ t1 + t2 , data = nuclearplants,
                            within = exactMatch(pr ~ pt, nuclearplants))

### Mahalanobis matching within subgroups, with a propensity score
### caliper:
fullmatch(match_on.examples$mh2 + caliper(match_on.examples$ps3, 1), data = nuclearplants)

### Alternative methods to matching without groups (exact matching)
m1 <- match_on(pr ~ t1 + t2, data=nuclearplants, within=exactMatch(pr ~ pt, nuclearplants))
m2 <- match_on(pr ~ t1 + t2 + strata(pt), data=nuclearplants)
# m1 and m2 are identical

m3 <- match_on(glm(pr ~ t1 + t2 + cost, data=nuclearplants,
                   family=binomial),
               data=nuclearplants,
               within=exactMatch(pr ~ pt, data=nuclearplants))
m4 <- match_on(glm(pr ~ t1 + t2 + cost + pt, data=nuclearplants,
                   family=binomial),
               data=nuclearplants,
               within=exactMatch(pr ~ pt, data=nuclearplants))
m5 <- match_on(glm(pr ~ t1 + t2 + cost + strata(pt), data=nuclearplants,
                   family=binomial), data=nuclearplants)
# Including `strata(foo)` inside a glm uses `foo` in the model as
# well, so here m4 and m5 are equivalent. m3 differs in that it does
# not include `pt` in the glm.
