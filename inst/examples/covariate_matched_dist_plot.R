data(nuclearplants)
### Full matching on a Mahalanobis distance.
fm1 <- fullmatch(pr ~ t1 + t2, data = nuclearplants)

library(ggplot2)
library(tidyverse)

nuclearplants$fm1 <- fm1

## The below description shows the following about differences within set on t2:
## half of the sets differ by less than 7
## 90% of the sets differ by less than 10.8
## The set with the largest difference had a difference of 13.
quantile(abs(t2_diffs$cov_diffs), c(0, .1, .25, .5, .9, 1))

## This next is a plot visualizing that distribution of within-stratum
## absolute differences

g0 <- covariate_matched_dist_plot(
  strata_indicator = fm1, trt = nuclearplants$pr,
  covariate = nuclearplants$t2
)
## Try the following to see the plot
## print(g0)
