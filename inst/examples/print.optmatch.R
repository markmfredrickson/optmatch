data(nuclearplants)
fm <- fullmatch(match_on(pr ~ t1 + t2, data = nuclearplants),
                data = nuclearplants)

print(fm)
print(fm, grouped = TRUE)
