options(repos = c(CRAN = "https://cran.rstudio.com/"))
devtools::release(args=c('--compact-vignettes=gs+qpdf') )
