em <- exactMatch(pr ~ pt, data=nuclearplants)
m1 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants)
stratumStructure(m1)
(subdims_em <- subdim(em))
m2 <- fullmatch(pr ~ t1 + t2, within=em, data=nuclearplants,
                mean.controls=pmin(1.5, subdims_em["controls",] / subdims_em["treatments",])
                )
stratumStructure(m2)
                
