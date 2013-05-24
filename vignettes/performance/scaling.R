library(optmatch)

Ks <- c(10, 20, 50, 100, 200, 300, 350, 400, 425, 450, 475, 500)
propensity <- rnorm(2 * max(Ks))
names(propensity) <- 1:(2 * max(Ks))

times <- lapply(Ks, function(k) {
  data <- propensity[1:(2 * k)]
  z <- rep(c(0,1), k)
  mdt <- system.time(mdist(z ~ data), gcFirst = TRUE) 
  mot <- system.time(match_on(x = data, z = z), gcFirst = TRUE)

  # for both runs, return the sum of the user and system time
  return(c(mdt[1] + mdt[2], mot[1] + mot[2]))
})

times <- do.call(rbind, times)
colnames(times) <- c("mdist", "match_on")
times <- data.frame(k = Ks, times)

save(file = "scaling.rda",
     times)

