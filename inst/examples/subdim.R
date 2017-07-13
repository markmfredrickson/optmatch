em <- exactMatch(pr ~ pt, data=nuclearplants)
m1 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants)
stratumStructure(m1)
(subdims_em <- subdim(em))
m2 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants,
                mean.controls=pmin(1.5, sapply(subdims_em, function(x) x[2]/x[1]))
                )
stratumStructure(m2)
                
