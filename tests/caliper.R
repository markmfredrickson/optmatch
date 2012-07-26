require("optmatch")

test <- function(t, m = "Error!") {
  if (!t) {
    stop(m)
  } 
}

shouldError <- function(expr, msg = "Exception should be thrown") {
  r <- try(expr, silent = T)
  if (!inherits(r, "try-error")) {
    stop(msg)  
  }
}

stripCall <- function(obj) {
  attr(obj, "call") <- NULL
  obj
}

### Basic test of caliper function ###
data(nuclearplants)

# use the Mahalanobis distance mdist method
result <- caliper(.2, pr ~ t1 + t2, data = nuclearplants)
test(inherits(result, "optmatch.dlist"), "Should be a optmatch object")

# checking that matrix is only Inf and 0
result <- caliper(.2, pr ~ t1 + t2, data = nuclearplants)
vr <- as.vector(result$m)
test(all(vr == Inf | vr == 0), "Not inf or 1")

### excluding some entries from having a caliper
# first, test edge condition of an error
#shouldError(caliper(.2, pr ~ t1 + t2, data = nuclearplants, exclude = c("XYZ")))

# test helper function
# result <- optmatch:::createWidths(
#   2,
#   mdist( pr ~ t1 + t2, data = nuclearplants)$m,
#   c("A" = 3))
#   
# test(result["A",] == 3)
# test(result["B",] == 2)

# test exclusion functionality
result <- caliper(.2, pr ~ t1 + t2, data = nuclearplants, exclude = c("A", "B"))
test(all(result$m["A",] == 0))
test(all(result$m["B",] == 0))

### Calipers should have same attributes as normal matches ###
a <-mdist(pr ~ t1 + t2, data = nuclearplants)
b <- caliper(.5, pr ~ t1 + t2, data = nuclearplants)
test(identical(attributes(stripCall(a)), attributes(stripCall(b))))
