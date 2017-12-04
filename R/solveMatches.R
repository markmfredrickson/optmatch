SolveMatches <- function(rownames, colnames, distspec, min.cpt,
                         max.cpt, tolerance, omit.fraction=NULL, matched.distances=FALSE)
{
  if (min.cpt <=0 | max.cpt<=0) {
    stop("inputs min.cpt, max.cpt must be positive")
  }

  if (!all(rownames %in% dimnames(distspec)[[1]])) {
    stop("input \'rownames\' may only contain row names of input \'distspec\'")
  }

  if (!all(colnames %in% dimnames(distspec)[[2]])) {
    stop("input \'rownames\' may only contain col. names of input \'distspec\'")
  }

  # distance must have a prepareMatching object
  if (!hasMethod("prepareMatching", class(distspec))) {
    stop("Argument \'distspec\' must have a \'prepareMatching\' method")
  }

  # convert the distspec into a cannonical matching specification with columns
  # treated, control, distance
  dm <- prepareMatching(distspec)

  rownames <- as.character(rownames)
  colnames <- as.character(colnames)

  rfeas <- length(unique(dm$treated))
  cfeas <- length(unique(dm$control))

  # If any controls were unmatchable, they were dropped by prepareMatching, and
  # positive `omit.fraction`'s need to be updated.
  if (cfeas < length(colnames) & is.numeric(omit.fraction) && omit.fraction >0) {
    original_number_to_omit <- omit.fraction*length(colnames)
    number_implicitly_omitted_already <- length(colnames) - cfeas
    omit.fraction <- (original_number_to_omit - number_implicitly_omitted_already)/cfeas
    # This can happen if the number to be omitted is less than the number of unmatchables
    if (omit.fraction <= 0) {
      omit.fraction <- NULL
    }
  }

  # ... and similarly in the case of negative `omit.fraction` if there were
  # treatments that couldn't be matched.
  if (rfeas < length(rownames) & is.numeric(omit.fraction) && omit.fraction <0) {
    original_number_to_omit <- -1*omit.fraction*length(rownames)
    number_implicitly_omitted_already <- length(rownames) - rfeas
    omit.fraction <- - (original_number_to_omit - number_implicitly_omitted_already)/rfeas
    # This can happen if the number to be omitted is less than the number of unmatchables
    if (omit.fraction >= 0) {
      omit.fraction <- NULL
    }
  }

  if (floor(min.cpt) > ceiling(max.cpt) | ceiling(1/min.cpt) < floor(1/max.cpt))
  {
    ans <- rep("NA",length(rownames)+length(colnames))
    names(ans) <- c(rownames, colnames)
    return(list(cells=ans, maxerr=NULL, distance=NULL))
  }

  # the next block of code, the dm <- ... is commented out as
  # dm is no longer a matrix. Completely unreachable entries may be a
  # problem later, but
  if (is.null(omit.fraction)) {
    f.ctls <- 1
    # dm <- matrix(dm[rfeas, cfeas], sum(rfeas), sum(cfeas),
    # dimnames=list(rownames[rfeas], colnames[cfeas]))
  } else {
    if (!is.numeric(omit.fraction) | omit.fraction <0 | omit.fraction > 1) {
      stop("omit.fraction must be null or between 0 and 1")
    }

    f.ctls <- 1-omit.fraction
    # dm <- matrix(dm[rfeas,], sum(rfeas), length(colnames),
    #              dimnames=list(rownames[rfeas], colnames))
  }


  if (any(rfeas) & any(cfeas))
  {
    old.o <- options(warn=-1)
    options(old.o)
    if(all(dm$distance == floor(dm$distance)) & any(dm$distance > 0)) #checking if all distances are integer
    {
      #integer version
      temp.with.nodes <- intSolve(dm, min.cpt, max.cpt, f.ctls) #note the new structure of temp.with.nodes

    }
    else #double precision
    {
      temp.with.nodes <- DoubleSolve(dm, rfeas, cfeas, min.cpt, max.cpt, tolerance, reso, f.ctls)
    }



  } else {
    temp <- 0 ; maxerr <- 0 ; dist <- 0
  }
  temp <- temp.with.nodes$temp
  temp$treated <- factor(temp$treated)
  temp$control <- factor(temp$control)
  ans <- rep(NA,length(rownames)+length(colnames))
  names(ans) <- c(rownames, colnames)

  matches <- solution2factor(temp)
  ans[names(matches)] <- matches

  #node.prices.d <- temp.with.nodes$node.price.ints / reso #check that reso = 1 when in integer case
  return(list(cells = ans, err = temp.with.nodes$maxerr, node.prices = temp.with.nodes$node.prices)) #could probably just add back translation to double solver?
}


DoubleSolve <- function(dm, rfeas, cfeas, min.cpt,
                        max.cpt, tolerance, reso, f.ctls)
{
  if (any(dm$distance > 0)) {
    reso <- (.Machine$integer.max/64 -2)/max(dm$distance)
  } else {
    reso <- min(.Machine$integer.max/64 -2, (sum(rfeas)+sum(cfeas))/tolerance)
  }

  if (tolerance>0 & sum(rfeas)>1 & sum(cfeas)>1) {
    reso <- min(reso, (sum(rfeas) + sum(cfeas) - 2)/tolerance)
  }

  #options(old.o) don't think this line is super important
  .matcher <- function(dm, toIntFunction, reso, min.cpt, max.cpt, f.ctls) {
    tmp <- dm
    tmp$distance <- toIntFunction(dm$distance * reso)
    # fmatch(tmp,
    #        max.row.units = ceiling(1/min.cpt),
    #        max.col.units = ceiling(max.cpt),
    #        min.col.units = max(1, floor(min.cpt)), f=f.ctls)

    obj <- intSolve(tmp, min.cpt, max.cpt, f.ctls)
    return(obj)
  }

  # fmatch returns a matrix with columns `treatment`, `control`, and `solution`
  # it also has a column `distance` with toIntFuction(dm * reso)

  temp.with.nodes <- .matcher(dm, floor, reso, min.cpt, max.cpt, f.ctls)
  #temp <- temp.with.nodes$temp

  if (any(is.na(temp.with.nodes$temp$solution))) {
    maxerr <- 0
  } else {
    maxerr <- sum(temp.with.nodes$temp$solution * dm$distance, na.rm = TRUE) -
      sum(temp.with.nodes$temp$solution * temp.with.nodes$temp$distance, na.rm = TRUE) / reso +
      (sum(rfeas) > 1 & sum(cfeas) > 1) *
      (sum(rfeas) + sum(cfeas) - 2 - sum(temp.with.nodes$temp$solution)) / reso
  }

  if (maxerr > tolerance)
  {
    temp1 <- temp.with.nodes
    temp2 <- .matcher(dm, round, reso, min.cpt, max.cpt, f.ctls)

    if  (sum(temp1$temp$solution * dm$distance, na.rm = TRUE) <= sum(temp2$temp$solution * dm$distance, na.rm = TRUE)) {
      temp.with.nodes <- temp1
    } else {
      temp.with.nodes <- temp2
    }

    maxerr <- sum(temp.with.nodes$temp$solution * dm$distance, na.rm = TRUE) -
      sum(temp1$temp$solution * temp.with.nodes$temp$distance, na.rm = TRUE)/reso +
      (max(1, sum(rfeas) - 1) + max(1, sum(cfeas) - 1) -
         (sum(rfeas) == 1 & sum(cfeas) == 1) - sum(temp1$temp$solution)) / reso
  }
  temp.with.nodes$node.prices <- temp.with.nodes$node.prices / reso
  temp.with.nodes$maxerr <- maxerr
  return(temp.with.nodes)

}

# a small helper function to turn a solution data.frame into a factor of matches
solution2factor <- function(s) {
  s2 <- s[s$solution == 1,]

  if (dim(s2)[1] == 0) {
    return(NULL)
  }

  # control units are labeled by the first treated unit to which they are connected
  # unlist(as.list(...)) was the best way I could find to make this into a vector, keeping names
  control.links <- unlist(as.list(by(s2, s2$control, function(x) { x[1,"treated"] })))

  # treated units are labeld by the label of the first control unit to which they are connected
  treated.links <- unlist(as.list(by(s2, s2$treated, function(x) { control.links[x[1, "control"]][1] })))

  # join the links
  return(c(treated.links, control.links))

}

intSolve <- function(dm, min.cpt, max.cpt, f.ctls)
{
  temp <- fmatch(dm, max.row.units = ceiling(1/min.cpt), max.col.units = ceiling(max.cpt), min.col.units = max(1, floor(min.cpt)), f=f.ctls)
  temp.extended <- temp

  indx <- temp.extended$control %in% temp.extended$control[which(temp.extended$treated == '(_Sink_)')] & temp.extended$treated == '(_End_)'
  sink.node.price.v <- temp.extended$reduced.cost[which(temp.extended$treated == '(_Sink_)')] - temp.extended$reduced.cost[indx]
  if(all(sink.node.price.v == sink.node.price.v[1]))
  {
    sink.node.price <- sink.node.price.v[1]
  }
  else
  {
    stop('unexpected mismatch of bookkeeping node prices')
  }

  temp <- temp[1:(dim(dm)[1]),]
  c(which(temp.extended$control == '(_End_)'), which(temp.extended$treated == '(_End_)'))

  node.prices.i <- c(-temp.extended$reduced.cost[c(which(temp.extended$control == '(_End_)'), which(temp.extended$treated == '(_End_)'))],0, sink.node.price)
  match.with.node.prices <- list()
  match.with.node.prices$temp <- temp
  match.with.node.prices$node.prices <- node.prices.i
  return(match.with.node.prices)
}

