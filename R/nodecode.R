#' @export
#' Assumes that new treatment and control nodes have already been found and categorized -- that is, of the 'new' nodes, we know whether they are treatment or control
#' Ignoring some of the potential larger implications for now...not entirely sure what the context of 'new' means here.
handle_new_nodes <- function(new.ts, new.cs, old.node.data)
{

}
# assuming no new nodes have been added
# want to return list of vectors node prices, named factor -- each vector should correspond to subproblem
# nt + nc + 2 for each subproblem
# can nodes appear in multiple subproblems?
# handles preparation and extraction of nodes
#' @export
prep_warm_nodes <- function(problems, old.node.data)
{
  .create.node.vecs <- function(problem)
  {
    browser()
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
    }
    # detect new nodes, just putting in some placeholder stuff for now
    # new.ts <- problem[!r%in% old.node.data$name,]
    # new.cs <- problem[,!c%in% old.node.data$name]
    # nodes.to.pass<- handle_new_nodes(new.ts, new.cs, old.node.data)
    return(nodes.to.pass)
  }

  node.list <- mapply(FUN = .create.node.vecs, problem = problems, SIMPLIFY = FALSE)
  return(node.list)
  #should return a list of vectors or data frames, with each entry being the warm start values for each subproblem

}

#builds node.data df from a solution
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
#' @export
assemble_node.data <- function(solutions)
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
  return(do.call(rbind, lapply(solutions, FUN = ff)))
}

#' @export
assemble_prob.data <- function(solutions, subproblemids, min.controls, max.controls, out.mean.controls, out.omit.fraction)
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
    new.df$mean.control[which(new.df$group == subproblemids)] <- out.mean.controls
    new.df$omit.fraction[which(new.df$group == subproblemids)] <- out.omit.fraction
  }

  return(new.df)
}
