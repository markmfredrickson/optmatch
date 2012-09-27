data(nuclearplants)
# as.matrix() used for convenient presentation only

### Caliper  of .2 pooled SDs in the propensity score
ppty <- glm(pr ~ . - (pr + cost), family = binomial(), data = nuclearplants)
ppty.dist <- match_on(ppty)

as.matrix(pptycaliper <- caliper(ppty.dist, width = .2)) 

### caliper on the Mahalanobis distance
as.matrix(caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3))

### Combining a Mahalanobis distance matching with a caliper
### of 1 pooled SD in the propensity score:
as.matrix(mhd.pptyc <- caliper(ppty.dist, width = 1) +
          match_on(pr ~ t1 + t2, data = nuclearplants))
pairmatch(mhd.pptyc)

### Excluding observations from caliper requirements:
as.matrix(caliper(match_on(pr ~ t1 + t2, data = nuclearplants), width = 3, exclude = c("A", "f")))

