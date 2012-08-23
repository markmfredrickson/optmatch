data(nuclearplants)
### Full matching on a Mahalanobis distance
mhd  <- mdist(pr ~ t1 + t2, data = nuclearplants)
( fm1 <- fullmatch(mhd) )
summary(fm1)
### Full matching with restrictions
( fm2 <- fullmatch(mhd, min=.5, max=4) )
summary(fm2)
### Full matching to half of available controls
( fm3 <- fullmatch(mhd, omit.fraction=.5) )
summary(fm3)
### Full matching within a propensity score caliper.
ppty <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
### Note that units without counterparts within the
### caliper are automatically dropped.
( fm4 <- fullmatch(mhd + caliper(mdist(ppty), 1)) )
summary(fm4)

### Propensity balance assessment. Requires RItools package.
if (require(RItools)) summary(fm4,ppty) 

### The order of the names in the match factor is unspecified.
### To ensure matching of matched set to obserservation, you
### must align it by observation name with your data set. 
cbind(nuclearplants, matches = fm4[row.names(nuclearplants)])
