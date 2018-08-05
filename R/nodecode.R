#' Function to handle preparation and extraction of node data. If information about nodes from a previously solved problem exists,
#' #' this function will extract that information and prepare it to be used in the following matching procedure. Currently is very much in development, not really ready for general use.
#' @export
prep_warm_nodes <- function(problems, old.node.data, old.prob.data)
{

  .create.node.vecs <- function(problem)
  {

    r <- rownames(problem)
    c <- colnames(problem)
    ns <- c(r,c)
    groups.of.nodes.to.pass <- old.node.data[old.node.data$name %in% ns, c("group")]
    gs <- unique(groups.of.nodes.to.pass)
    if(length(gs) > 1)
    {
      warning("groups/subproblems have changed")
    }
    else
    {
      indx <- old.node.data$name %in% ns | (old.node.data$name == "(_End_)" & old.node.data$group == gs[1] | old.node.data$name == "(_Sink_)" & old.node.data$group == gs[1])
      nodes.to.pass <- old.node.data[indx, c("price")]
      names(nodes.to.pass) <- old.node.data[indx, "name"]
      reso.m <- old.prob.data[old.prob.data$group == gs, "reso"]
    }
    ###***** start of code that is not quite done: handling new nodes that have been introduced into a problem
    #consider case where the nodes in a subproblem might all be new, won't worry about it for now
    new.ts <- rownames(problem)[!rownames(problem) %in% old.node.data$name]
    new.cs <- colnames(problem)[!colnames(problem) %in% old.node.data$name]
    if(length(new.ts) > 0 | length(new.cs) > 0)
    {
      prices.new <- numeric(length = (length(r) + length(c) + 2))
      names(prices.new) <- c(r, c, "(_End_)", "(_Sink_)")
      prices.new[names(nodes.to.pass)] <- nodes.to.pass
      price.suggestions.c <- NULL
      price.suggestions.t <- NULL

      .handle.new.c <- function(new.c)
      {

        distvec <- as.vector(problem[,new.c]) #should return a subsetted df with just the info we want to examine
        prices <- ifelse(is.na(nodes.to.pass[r]), NA, nodes.to.pass[r])
        if(!is.na(reso.m))
        {
          distvec <- round(reso.m * distvec)
        }
        new.control <- function(dist.val, node.price)
        {

          if(is.infinite(dist.val))
          {
            #terminate and do something
            return(-Inf)
          }
          else if(is.na(node.price))
          {
            return(NA)
            #also terminate and do something
          }
          else
          {
            comp.val <- node.price - dist.val
            return(comp.val)
            #handle updating other metadata elsewhere
          }
        }
        A <- mapply(new.control, dist.val = distvec, node.price = prices)
        maxa <- max(A) #update the max value as other values are calculated? rather than having to go find a max
        return(maxa - 1)
      }

      old.nas <- length(new.cs)
      if(length(new.cs))
      {
        price.suggestions.c <- sapply(new.cs, .handle.new.c)
        if(!is.na(reso.m))
        {
          price.suggestions.c <- price.suggestions.c / reso.m
        }
        prices.new[names(price.suggestions.c)] <- price.suggestions.c

      }

      .handle.new.t <- function(new.t)
      {
        distvec <- as.vector(problem[new.t,])
        prices <- ifelse(is.na(nodes.to.pass[c]), NA, nodes.to.pass[c])
        #distvec should be some multiple of length(prices)
        if(!is.na(reso.m))
        {
          distvec <- round(reso.m * distvec)
        }

        new.treatment <- function(dist.val, node.price)
        {
          if(is.infinite(dist.val))
          {
            #terminate and do something
            return(Inf)
          }
          else if(is.na(node.price))
          {
            return(NA)
            #also terminate and do something
          }
          else
          {
            comp.val <- dist.val + node.price
            return(comp.val)

          }
        }
        A <- mapply(new.control, dist.val = distvec, node.price = prices)
        mina <- min(A) #update the max value as other values are calculated? rather than having to go find a max
        return(mina + 1)
      }
      old.nas.t <- length(new.ts)
      if(length(new.ts))
      {
        price.suggestions.t <- sapply(new.ts, .handle.new.t)
        if(!is.na(reso.m))
        {
          price.suggestions.t <- price.suggestions.t / reso.m
        }
        prices.new[names(price.suggestions.t)] <- price.suggestions.t
      }

      if(!is.null(price.suggestions.c) & (sum(is.na(price.suggestions.c)) < old.nas.c & sum(is.na(price.suggestions.c)) > 0))
      {
        new.cs <- new.cs[!is.na(price.suggestions.c)]
        price.suggestions.c <- sapply(new.cs, .handle.new.c)
        prices.new[names(price.suggestions.c)] <- price.suggestions.c
      }
      if(!is.null(price.suggestions.t) & (sum(is.na(price.suggestions.t)) < old.nas.t & sum(is.na(price.suggestions.t)) > 0))
      {
        new.ts <- new.ts[!is.na(price.suggestions.t)]
        price.suggestions.t <- sapply(new.ts, .handle.new.t)
        prices.new[names(price.suggestions.t)] <- price.suggestions.t
      }
      nodes.to.pass <- ifelse(is.infinite(prices.new) | is.na(prices.new), 0, prices.new)
    } ### *** end of 'new node handling code'

    return(nodes.to.pass)
  }

  node.list <- mapply(FUN = .create.node.vecs, problem = problems, SIMPLIFY = FALSE)
  return(node.list)
  #should return a list of vectors or data frames, with each entry being the warm start values for each subproblem

}

#builds node.data (as defined in Optmatch class definition) data frame from a subproblem matching solution
#' @export
build_node_data <- function(temp.extended, subproblemid, treatment.names, control.names)
{


  nnodes <- c(as.character(temp.extended$treated[which(temp.extended$control == '(_End_)')]), as.character(temp.extended$control[which(temp.extended$treated == '(_End_)')]), "(_End_)", "(_Sink_)")
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
  node.prices.i <- c(-temp.extended$reduced.cost[c(which(temp.extended$control == '(_End_)'), which(temp.extended$treated == '(_End_)'))],0, sink.node.price)

  node.data <- data.frame(name = nnodes, price = node.prices.i)

  .determine_group <- function(node.name)
  {
      if(node.name == "(_End_)" | node.name == "(_Sink_)")
      {
        return(NA) #bookkeeping
      }
      if(all(node.name %in% temp.extended$treated))
      {
        return(as.logical(TRUE)) #treatment
      }
      else
      {
        return(as.logical(FALSE)) #control
      }
  }
  node.data$contrast.group <- sapply(nnodes, FUN = .determine_group)
  # node.data <- data.frame(name = c(nnodes, msg.n), price = c(node.prices.i, msg.p), contrast.group = c(cg, msg.c), group <- c(rep(subproblemid, length(nnodes)), msg.s)
  node.data$group <- subproblemid
  return(node.data)
}

# Combines all information from node.data data frames across all subproblems into a single node.data frame for the entire original problem
#' @export
assemble_node.data <- function(solutions, treated = NULL, control = NULL)
{
  ff <- function(x){
    if(is.null(x$node.data))
    {
      return(data.frame())
    }
    else
    {
      return(x$node.data)
    }
  }
  df <- do.call(rbind, lapply(solutions, FUN = ff))
  #browser()
  if(!is.null(treated) & !is.null(control))
  {
    dropped.t <- treated[!treated %in% df$name]
    dropped.c <- control[!control %in% df$name]
    #df2 <- data.frame(name = c(dropped.t, dropped.c), price = NA, contrast.group = c(rep(TRUE, rep(length(dropped.t))), rep(FALSE, length(dropped.c))), group = NA)
    df2 <- data.frame(name = c(dropped.t, dropped.c), price = rep(NA, length(c(dropped.t, dropped.c))))
    df2$contrast.group <- c(rep(TRUE, rep(length(dropped.t))), rep(FALSE, length(dropped.c)))
    df2$group <- rep(NA, length(c(dropped.t, dropped.c)))
    attr(df, "dropped.nodes") <- df2
  }

  rownames(df) <-NULL
  return(df)
}

#compiles all prob.data data frames across all subproblems into a single prob.data data frame for the entire original matching problem
#' @export
assemble_prob.data <- function(solutions, subproblemids = NA, min.controls = NA, max.controls = NA, out.mean.controls = NA, out.omit.fraction = NA)
{
  ff <- function(x){
    if(is.null(x$prob.data))
    {
      return(data.frame())
    }
    else
    {
      return(x$prob.data)
    }
  }
  new.df <- do.call(rbind, lapply(solutions, FUN = ff))

  if(nrow(new.df))
  {
    new.df$min.control[which(new.df$group == subproblemids)] <- min.controls
    new.df$max.control[which(new.df$group == subproblemids)] <- max.controls
    if(is.null(new.df$mean.control))
    {
      new.df$mean.control <- NA
    }
    new.df$mean.control[which(new.df$group == subproblemids)] <- out.mean.controls
    new.df$omit.fraction[which(new.df$group == subproblemids)] <- out.omit.fraction
  }

  return(new.df)
}
