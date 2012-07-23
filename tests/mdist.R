require("optmatch")

test <- function(t, m = "Test Failed!") {
  if (!t) {
    stop(m)
  } 
}

### Data ###

data("nuclearplants")

fmla <- pr ~ t1 + t2 + pt
test.glm <- glm(fmla, family=binomial(), data=nuclearplants)

### Basic Tests ###

result.glm <- mdist(test.glm)
result.glm2 <- mdist(test.glm, ~pt)

test(inherits(result.glm, "optmatch.dlist"), "Should be a optmatch object")
test(length(result.glm) == 1)
test(length(result.glm2) == 2)

result.optmatch.dlist <- mdist(result.glm)
test(inherits(result.optmatch.dlist, "optmatch.dlist"), "Should be a optmatch object")
test(length(result.optmatch.dlist) == 1)

result.fmla <- mdist(fmla, data = nuclearplants)
test(inherits(result.fmla, "optmatch.dlist"), "Should be a optmatch object")
test(length(result.fmla) == 1)

result.fmla2 <- mdist(fmla, ~pt, data = nuclearplants)
test(inherits(result.fmla2, "optmatch.dlist"), "Should be a optmatch object")
test(length(result.fmla2) == 2)

### Function Tests ###

# first, a simpler version of scalar diffs
sdiffs <- function(treatments, controls) {
  abs(outer(treatments$t1, controls$t1, `-`))
}

result.function <- mdist(sdiffs, pr ~ 1, nuclearplants)
test(all(dim(result.function) == c(10,22)), "Function not returning right sized data")

test(identical(optmatch:::parseFmla(y ~ a | group), lapply(c("y", "a", "group"), as.name)))
test(identical(optmatch:::parseFmla(y ~ a), c(as.name("y"), as.name("a"), NULL)))

test(!is.null(rownames(result.function$m)) && all(rownames(result.function$m) %in% rownames(nuclearplants[nuclearplants$pr == 1,])))

result.function.a <- mdist(sdiffs, pr ~ 1 | pt, nuclearplants)
result.function.b <- mdist(sdiffs, pr ~ pt, nuclearplants)

test(identical(result.function.a, result.function.b), "Two ways to specify groupings")

shouldError <- function(expr, msg = "Exception should be thrown") {
  r <- try(expr, silent = T)
  if (!inherits(r, "try-error")) {
    stop(msg)  
  }
}

shouldError(mdist(sdiffs, pr ~ pt + t1, nuclearplants))

# the fun part, making a dlist when there are multiple groups

test(length(mdist(sdiffs, pr ~ pt, nuclearplants)) == 2, "Not enough groups")

### Using mad() instead of sd() for GLM distances

result <- mdist(glm(pr ~ t1 + t2 + cost, data = nuclearplants, family = binomial()))

# this is an odd test, but a simple way to make sure mad is running, not SD(). 
# I would like a better test of the actual values, but it works
test(mean(result$m) > 2)


### mdist() should informatively complain if passed a numeric vector.
### (This may change in the future.)

shouldError(mdist(test.glm$linear.predictor))


### Stratifying by a pipe (|) character in formulas

main.fmla <- pr ~ t1 + t2
strat.fmla <- ~ pt
combined.fmla <- pr ~ t1 + t2 | pt

result.main <- mdist(main.fmla, structure.fmla = strat.fmla, data = nuclearplants)
result.combined <- mdist(combined.fmla, data = nuclearplants)

test(identical(result.main, result.combined))

### Informatively insist that one of formulas specify the treatment group
shouldError(mdist(~t1+t2, structure.fmla=~pt, data=nuclearplants))
test(identical(mdist(pr~t1+t2, structure.fmla=~pt, data=nuclearplants),
               mdist(~t1+t2, structure.fmla=pr~pt, data=nuclearplants))
     )
### Finding "data" when it isn't given as an argument
### Caveats:
### * data's row.names get lost when you don't pass data as explicit argument;
### thus testing with 'all.equal(unlist(<...>),unlist(<...>))' rather than 'identical(<...>,<...>)'.
### * with(nuclearplants, mdist(fmla)) bombs for obscure scoping-related reasons,
### namely that the environment of fmla is the globalenv rather than that created by 'with'.
### This despite the facts that identical(fmla,pr ~ t1 + t2 + pt) is TRUE and that
### with(nuclearplants, mdist(pr ~ t1 + t2 + pt)) runs fine.
### But then with(nuclearplants, lm(fmla)) bombs too, for same reason, so don't worry be happy.
attach(nuclearplants)
test(all.equal(unlist(result.fmla),unlist(mdist(fmla))))
test(all.equal(unlist(result.main),unlist(mdist(main.fmla, structure.fmla=strat.fmla))))
test(all.equal(unlist(result.combined),unlist(mdist(combined.fmla))) )
detach("nuclearplants")
test(identical(fmla,pr ~ t1 + t2 + pt))
test(all.equal(unlist(result.fmla),unlist(with(nuclearplants, mdist(pr ~ t1 + t2 + pt)))))
test(identical(combined.fmla, pr ~ t1 + t2 | pt))
test(all.equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2 | pt)))))
test(all.equal(unlist(result.fmla), unlist(with(nuclearplants[-which(names(nuclearplants)=="pt")],
                                           mdist(update(pr ~ t1 + t2 + pt,.~.-pt + nuclearplants$pt))
                                           )
                                      )
          )
     )
test(all.equal(unlist(result.combined), unlist(with(nuclearplants, mdist(pr ~ t1 + t2, structure.fmla=strat.fmla)))))

### bigglm method
if (require('biglm')) {
bgps <- bigglm(fmla, data=nuclearplants, family=binomial() )
shouldError(mdist(bgps, structure.fmla=pr ~ 1))
shouldError(mdist(bgps, data=nuclearplants))
result.bigglm1 <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants)
result.bigglm1a <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
                        standardization.scale=sd)
result.bigglm1b <- mdist(bgps, structure.fmla=pr ~ 1, data=nuclearplants,
                        standardization.scale=NULL)
result.bigglm2 <- mdist(bgps, structure.fmla=pr ~ pt, data=nuclearplants)
}

### Jake found a bug 2010-06-14
### Issue appears to be a missing row.names/class

absdist1 <- mdist(sdiffs, structure.fmla = pr ~ 1|pt, data = nuclearplants)
test(length(pairmatch(absdist1)) > 0)

### Check that distances combine as they should
### (a joint test of mdist and Ops.optmatch.dlist)
### Distances without subclasses:
test(inherits(result.glm + result.fmla, "optmatch.dlist"),
     "Should be a optmatch object")
test(inherits(result.glm + result.function, "optmatch.dlist"),
     "Should be a optmatch object")
if (require("biglm"))
test(inherits(result.glm + result.bigglm1, "optmatch.dlist"),
     "Should be a optmatch object")


### Distances embodying subclassification:
test(inherits(result.glm2 + result.fmla2, "optmatch.dlist"),
     "Should be a optmatch object")

test(inherits(result.glm2 + result.function.a, "optmatch.dlist"),
     "Should be a optmatch object")

tmp <- result.glm2[1]
class(tmp) <- c("optmatch.dlist", 'list')

shouldError(inherits(tmp + result.function.a, "optmatch.dlist"),
     "Should be a optmatch object")

if (require("biglm"))
test(inherits(result.glm2 + result.bigglm2, "optmatch.dlist"),
     "Should be a optmatch object")
