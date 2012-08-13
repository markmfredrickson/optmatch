fullmatch <- function(distance,  
    min.controls = 0, 
    max.controls = Inf, 
    omit.fraction = NULL, 
    tol = .001) {

  ### Checking Input ###
  
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification object")      
  }

  # we expect the following to be defined for the distance object
  # put any functions in this list that are called directly on distance
  methods <- c("dim", "dimnames", "prepareMatching", "subproblems",
    "is.numeric")
  dist.class <- class(distance)
  lapply(methods, function(m) {
    if (!hasMethod(m, dist.class)) {
      # skip the FUN = ... in the call stack
      stop(paste("argument \'distance\' must have a", m, "method."), call. = F)  
    }
  })

  dnms <- dimnames(distance)
  if (is.null(dnms) | is.null(dnms[[1]]) | is.null(dnms[[2]])) {
    stop("argument \'distance\' must have dimnames") 
  }
  
  if (any(duplicated(unlist(dnms)))){ 
    stop("dimnames of argument \'distance\' contain duplicates") 
  }

  nmtrt <- dnms[[1]]
  nmctl <- dnms[[2]]

  # note: this next _should_ be unnecessary, the objects should do this
  # but better safe than sorry
  if (!identical(dim(distance), c(length(nmtrt), length(nmctl)))) {
    stop("argument \'distance\' dimensions do not match row and column names")  
  }

  # distances and other args must be numeric
  if (!is.numeric(distance)) {
    stop("argument \'distance\' must be numeric")  
  }
  if (!is.numeric(min.controls)) {
    stop("argument \'min.controls\' must be numeric")  
  }
  if (!is.numeric(max.controls)) {
    stop("argument \'max.controls\' must be numeric")  
  }
  if (!is.null(omit.fraction)) {
    if (any(abs(omit.fraction) > 1, na.rm = TRUE) | !is.numeric(omit.fraction)) { 
      stop("omit.fraction must be NULL or numeric between -1 and 1") 
    }
  }


  # problems is guaranteed to be a list of DistanceSpecifictions
  # it may only have 1 entry
  problems <- findSubproblems(distance)

  # the number of problems should match the argument lengths for
  # min, max, and omit

  np <- length(problems)
  if (length(min.controls) > 1 & np != length(min.controls)) {
    stop(paste("Length of \'min.controls\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  if (length(max.controls) > 1 & np != length(max.controls)) {
    stop(paste("Length of \'max.controls\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  if (!is.null(omit.fraction) & length(omit.fraction) > 1 & np !=
    length(omit.fraction)) {
    stop(paste("Length of \'omit.fraction\' arg must be same ",
              "as number of subproblems [", np, "]", sep = ""))  
  }
  
  # reset the arguments to be the right length if they are not
  if (length(min.controls) == 1) {
    min.controls <- rep(min.controls, np)   
  }
  if (length(max.controls) == 1) {
    max.controls <- rep(max.controls, np)   
  }

  if (is.null(omit.fraction)) {
    omit.fraction <- NA  
  }
  if (length(omit.fraction) == 1) {
    omit.fraction <- rep(omit.fraction, np)   
  }

  if (any(omit.fraction > 0 & max.controls <= .5, na.rm=TRUE)) {
      stop("positive \'omit.fraction\' with \'max.controls\' <= 1/2 not permitted")    
  }

  if (any(omit.fraction < 0 & min.controls >= 2, na.rm=TRUE)) {
      stop("negative \'omit.fraction\' with \'min.controls\' >= 2 not permitted")
  }
  
  total.n <- sum(dim(distance))

  TOL <- tol * total.n

  # a helper to handle a single matching problem. all args required. 
  # input error checking happens in the public fullmatch function.
  .fullmatch <- function(d, mnctl, mxctl, omf) {
    
    # if the subproblem is completely empty, short circuit
    if (length(d) == 0 || all(is.infinite(d))) {
      x <- dim(d)
      cells.a <- rep(NA, x[1])
      cells.b <- rep(NA, x[2])
      names(cells.a) <- rownames(d)
      names(cells.b) <- colnames(d)
      tmp <- list(cells = c(cells.a, cells.b), maxerr = -1)
      return(tmp) 
    }

    ncol <- dim(d)[2]
    nrow <- dim(d)[1]

    tol.frac <- (nrow + ncol - 2)/(total.n - 2 * np)
  
    # if omf is specified (i.e. not NA), see if is greater than 0
    # if omf is not specified, check to see if mxctl is > .5
    if (switch(1 + is.na(omf), omf > 0,  mxctl > .5)) {
      maxc <- min(mxctl, ncol)
      minc <- max(mnctl, 1/nrow)
      omf.calc <- omf
  
    } else {
      maxc <- min(1/mnctl, ncol)
      minc <- max(1/mxctl, 1/nrow)
      omf.calc <- -1 * omf
    }

    temp <- SubDivStrat(rownames = rownames(d),
                        colnames = colnames(d),
                        distspec = d, 
                        max.cpt = maxc,
                        min.cpt = minc,
                        tolerance = TOL * tol.frac,
                        omit.fraction = if(!is.na(omf)) { omf.calc }) # passes NULL for NA

    return(temp)
  }
  
  solutions <- mapply(.fullmatch, problems, min.controls, max.controls, omit.fraction, SIMPLIFY = FALSE)

  return(makeOptmatch(distance, solutions, match.call()))
}


