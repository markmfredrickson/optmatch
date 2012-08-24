data(nuclearplants) 

### Pair matching on a Mahalanobis distance 
mhd <- mdist(pr ~ t1 + t2, data = nuclearplants) 
( pm1 <- pairmatch(mhd) ) 
summary(pm1) 

### Pair matching within a propensity score caliper.  
ppty <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants) 
( pm2 <- pairmatch(mhd + caliper(mdist(ppty), 2)) ) 
summary(pm2)

### Propensity balance assessment. Requires RItools package.
if(require(RItools)) summary(pm2, ppty)

### 1:2 matched triples
tm <- pairmatch(mhd, controls = 2)
summary(tm)

### Creating a data frame with the matched sets attached.
### mdist(), caliper() and the like cooperate with pairmatch()
### to make sure observations are in the proper order:
all.equal(names(tm), row.names(nuclearplants))
### So our data frame including the matched sets is just
cbind(nuclearplants, matches=tm)

### In contrast, if your matching distance is an ordinary matrix
### (as earlier versions of optmatch required), you'll
### have to align it by observation name with your data set. 
cbind(nuclearplants, matches = tm[row.names(nuclearplants)])

