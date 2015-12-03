data(nuclearplants)
fm <- fullmatch(pr ~ t1 + t2, data = nuclearplants)

print(fm)
print(fm, grouped = TRUE)
