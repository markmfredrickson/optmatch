require("optmatch")
library("splines")
data(nuclearplants)

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

### Basic Tests
# Takes and returns a data frame
test(is.data.frame(fill.NAs(data.frame(1))))

# A formula alone is not allowed
shouldError(fill.NAs(y ~ x))

# takes a formula and a data.frame, returns a data frame
result <- fill.NAs(pr ~ cost, nuclearplants)
test(is.data.frame(fill.NAs(pr ~ cost, nuclearplants))) # no missingness

# simple calls should be equivalent to model.frame
test(length(result) == 2)

# Adds additional columns for missing data indicators
sample.df <- data.frame(a = 1:100, b = 100:1, c = rep(c(1,2, NA, 3, 4), 20))
result <- fill.NAs(sample.df)
test(length(colnames(result)) == 4)

# the last column should be TRUE every 3 unit
test(identical(result[[4]], rep(c(F, F, T, F, F), 20)))

# column name should be c.NA
test(identical(colnames(result)[4], "c.NA"))

# for variables encapsulated in functions, only the variable should be expanded into a NA column
nuclear.missing <- nuclearplants
nuclear.missing$cost[sample(length(nuclear.missing$cost), 5)] <- NA
result <- fill.NAs(pr ~ ns(cost, df = 3), nuclear.missing)
test(length(result) == 5, "Expanding only the cost column, not the splines") # pr, 3 cost splines, cost.NA
test(colnames(result)[1] == "pr", "response should be first column")

### Imputation
# imputes values for main columns, not splines
# test(mean(result[[2]][,3][!result$cost.NA] != result[[2]][1,3]))

#### Using model.matrix under the hood makes every column numeric.
# logical imputation
#logical.df <- data.frame(y = c(1,2,3,4,5,6,7,8), l = c(T,T,T,F,F, NA, T,F))
#result.logical <- fill.NAs(logical.df)

#test(class(result.logical$l) == "logical")
#test(!(any(is.na(result.logical))))
#test(result.logical[[1]][6] == TRUE)

# factor imputation
#factor.df <- data.frame(f = factor(c("a", "b", "b", NA, "c")))
#result.factor <- fill.NAs(factor.df)

#test(class(result.factor[[1]]) == "factor")
#test(!(any(is.na(result.factor))))
#test(result.factor[[1]][4] == ".NA")

## ordered factor imputation
#ordered.numbers <- c(1,3,2,2,3,NA,2,4,1,1,4)
#ordered.df <- data.frame(o = ordered(ordered.numbers))
#result.ordered <- fill.NAs(ordered.df)

#test(inherits(result.ordered[[1]], "ordered"))
#test(!(any(is.na(result.ordered))))
#test(result.ordered[[1]][6]== 2)

## Matrices aren't data frames, but they should be valid
data(nuclearplants)
h <- dim(nuclearplants)[1]
w <- dim(nuclearplants)[2]
np.missing <- t(apply(
  cbind(nuclearplants, missing = rbinom(h, prob = .1, size = w)), 1, 
  function(row) { row[sample(w, row["missing"])] <- NA; row })) 

test(is.data.frame(fill.NAs(pr ~ t1 + t2, data = np.missing)))
test(is.data.frame(fill.NAs(np.missing)))

## Result should be directly passable to lm()
n2 <- nuclearplants
n2[1,"t1"] <- NA

imputed.fmla <- fill.NAs(cost ~ log(t1), data = n2)
imputed.frame <- fill.NAs(n2)
m1 <- lm(imputed.fmla)
m2 <- lm(cost ~ log(t1) + t1.NA, data = imputed.frame)

# for some reason log(t1) appears as `log(t1)`. I strip these
# out and treat the results as equal otherwise
test(identical(gsub("`", "", names(m1$coef)), names(m2$coef)))

## right number of columns if 2 of the same variable used
imputed.fmla <- fill.NAs(cost ~ log(t1) + sqrt(t1), data = n2)
test(dim(imputed.fmla)[2] == 4)

#### Do not impute response, only covariates
naresponse.df <- data.frame(Y = c(1, 2, 3, NA, 5), X = c(10, 20, NA, 40, 50))
imputed.response <- fill.NAs(Y ~ X, naresponse.df)
test(any(is.na(imputed.response$Y)))
test(!any(is.na(imputed.response$X)))

#### Transform, then impute ####
#### turning off tests for now. the strategy is to use model.matrix before
#### imputing
transform.df <- data.frame(Y = c(1,2,3,4,5), X1 = c(2,2,4, NA, 4), X2 = c(NA, 10, 20, 30, NA))
imputed.transform <- fill.NAs(Y ~ X1 * X2, data = transform.df)
# should have 6 columns Y, X1, X2, X2:X3, X1.NA, and X2.NA
test(dim(imputed.transform)[2] == 6) 
test(identical(imputed.transform$X1 , c(2,2,4,3,4)))
test(identical(imputed.transform$X2 , c(20, 10, 20, 30, 20)))
test(all(imputed.transform["X1:X2"] == c(50, 20, 80, 50, 50)))

i2.transform <- fill.NAs(Y ~ X1, data = transform.df)
test(length(i2.transform) == 3, "bug with adding .na columns for things not in the formula")

