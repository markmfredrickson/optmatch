library(optmatch)

exdir <- system.file("examples", package = "optmatch")

exs <- Sys.glob(file.path(exdir, "*"))

for (e in exs) {
  # jump through a few extra hoops to ensure examples don't conflict on vars, etc.
  etmp <- new.env()
  eval(envir = etmp,
       source(e, echo = T, local = TRUE, max.deparse.length = -1))
}
