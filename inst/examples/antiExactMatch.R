data(nuclearplants)

# force entries to be within the same factor:
em <- fullmatch(exactMatch(pr ~ pt, data = nuclearplants), data = nuclearplants)
table(nuclearplants$pt, em)

# force treated and control units to have different values of `pt`:
z <- nuclearplants$pr
names(z) <- rownames(nuclearplants)
aem <- fullmatch(antiExactMatch(nuclearplants$pt, z), data = nuclearplants)
table(nuclearplants$pt, aem)
