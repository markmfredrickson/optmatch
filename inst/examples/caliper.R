data(nuclearplants)


### Caliper of 100 MWe on plant capacity
caliper(match_on(pr~cap, data=nuclearplants, method="euclidean"), width=100)

### Caliper of 1/2 a pooled SD of plant capacity
caliper(match_on(pr~cap, data=nuclearplants), width=.5)

### Caliper  of .2 pooled SDs in the propensity score
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
ppty.dist <- match_on(ppty)

pptycaliper <- caliper(ppty.dist, width = .2)

### caliper on the Mahalanobis distance
caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3)

### Combining a Mahalanobis distance matching with a caliper
### of 1 pooled SD in the propensity score:
mhd.pptyc <- caliper(ppty.dist, width = 1) +
          match_on(pr ~ t1 + t2, data = nuclearplants)
pairmatch(mhd.pptyc, data = nuclearplants)

### Excluding observations from caliper requirements:
caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3, exclude = c("A", "f"))
