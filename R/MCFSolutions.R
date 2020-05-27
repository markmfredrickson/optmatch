setOldClass("data.frame")
####################################################################
########  Classes for storing information about solutions of #######
########  Min-Cost-Flow representations of matching problems #######
#######  (See vignette "MCFSolutions" for class descriptions.) #####
####################################################################
setClass("SubProbInfo", contains="data.frame",
         prototype=
             prototype(data.frame(groups=character(1), flipped=NA, hashed_dist=NA_character_,
                              resolution=1, lagrangian_value=NA_real_,
                              dual_value=NA_real_, feasible=NA, exceedance=NA_real_, stringsAsFactors=FALSE)
                  )
         )
setValidity("SubProbInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:8]==
             c("groups","flipped", "hashed_dist","resolution","lagrangian_value","dual_value", "feasible", "exceedance")))
        errors  <- c(errors,
                     'Cols 1-8 should be:\n\t c("groups","flipped", "hashed_dist","resolution","lagrangian_value","dual_value", "feasible", "exceedance")')
    if (!all(vapply(object[c(1,3)], is.character, logical(1))))
        errors  <- c(errors,
                     'Cols 1,3 should have type character.')
    if (!all(vapply(object[c(4:6,8)], is.double, logical(1))))
        errors  <- c(errors,
                     'Cols 4-6, 8 should have type double.')
    if (!all(vapply(object[c(2,7)], is.logical, logical(1))))
        errors  <- c(errors,
                     'Col 2,7 should have type logical.')
    if (anyDuplicated(object[['groups']]))
        errors  <- c(errors,
                     'Duplicates in "groups", or subproblems with same name.')
    if (length(errors)==0) TRUE else errors  
})
setClass("NodeInfo", contains="data.frame",
         prototype=
             prototype(data.frame(name = character(0), price=double(0),
                                  upstream_not_down=logical(0), supply=integer(0),
                                  groups=factor(), stringsAsFactors=FALSE)
                       )
         )
setValidity("NodeInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:5]==
             c("name", "price", "upstream_not_down", "supply", "groups")))
        errors  <- c(errors,
                     'Cols 1-5 should be:\n\t c("name", "price", "upstream_not_down", "supply", "groups")')
    if (!is.character(object[['name']]))
        errors  <- c(errors,
                     'Col "name" should have type character.')
    if (!is.factor(object[['groups']]))
        errors  <- c(errors,
                     'Col "groups" should have type factor.')
    if (!is.numeric(object[['price']]))
        errors  <- c(errors,
                     'Col "price" should be a numeric.')
    if (!is.logical(object[['upstream_not_down']]))
        errors  <- c(errors,
                     'Col "upstream_not_down" should have type logical.')
    if (!is.integer(object[['supply']]))
        errors  <- c(errors,
                     'Col "supply" should have type integer.')    
    if (length(errors)==0) TRUE else errors      
})

setClass("ArcInfo", slots=c(matches="data.frame", bookkeeping="data.frame"),
         prototype=
             prototype(matches=data.frame(groups = factor(), upstream = factor(), 
                                          downstream = factor(), stringsAsFactors=FALSE),
                       bookkeeping=data.frame(groups = factor(), start = factor(), 
                                              end = factor(), flow=integer(0),
                                              capacity=integer(0), stringsAsFactors=FALSE)
                       )
         )
setValidity("ArcInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object@matches)[1:3]==
             c("groups", "upstream",  "downstream")))
        errors  <- c(errors,
                     '@matches cols 1-3 should be:\n\t c("groups", "upstream",  "downstream")')
    if (!all(vapply(object@matches, is.factor, logical(1))==TRUE))
        errors  <- c(errors,
                     'All columns of @matches should have type factor.')
    if (!all(colnames(object@bookkeeping)[1:5]==
             c("groups", "start",  "end",  "flow", "capacity")))
        errors  <- c(errors,
                     '@bookkeeping cols 1-5 should be:\n\t c("groups", "start",  "end",  "flow", "capacity")')
    if (!all(vapply(object@bookkeeping[1:3], is.factor, logical(1))))
        errors  <- c(errors,
                     '@bookkeeping cols 1-3 should have type factor.')
    if (!all(vapply(object@bookkeeping[c('flow', "capacity")], is.integer, logical(1))))
        errors  <- c(errors,
                     '@bookkeeping cols "flow", "capacity" should have type integer.')
    if (!all(object@bookkeeping[['flow']]>=0L))
        errors  <- c(errors,
                     '@bookkeeping "flow" values should be nonnegative.')
    if (!all(object@bookkeeping[['flow']]<=object@bookkeeping[['capacity']]))
        errors  <- c(errors,
                     'in @bookkeeping, flow can be no greater than capacity.')
    if (length(errors)==0) TRUE else errors      
})


setClass("MCFSolutions", slots=c(subproblems='SubProbInfo',nodes='NodeInfo',
                                 arcs='ArcInfo'),
         prototype = prototype(subproblems=new('SubProbInfo'), nodes=new('NodeInfo'),
                               arcs=new('ArcInfo'))
         )
setValidity("MCFSolutions", function(object){
    errors  <- character(0)
    ## Each of factors object@arcs@matches$upstream, object@arcs@matches$downstream,
    ## object@arcs@bookkeeping$start and object@arcs@bookkeeping$end must have the same
    ## levels set, namely node.labels(object) (i.e. row.names(object@nodes) ). 
    if (length(xtralevs  <- setdiff(levels(object@arcs@matches[['upstream']]),
                                    node.labels(object)
                                    )
               )
        )
        errors  <- c(errors,
                     paste("Arcs' upstream nodes not listed in nodes table, e.g.",
                           paste(head(xtralevs,2), collapse=", "), "."
                           )
                     )
    if (length(xtralevs  <- setdiff(levels(object@arcs@matches[['downstream']]),
                                    node.labels(object)
                                    )
               )
        )
        errors  <- c(errors,
                     paste("Arcs' downstream nodes not listed in nodes table, e.g.",
                           paste(head(xtralevs,2), collapse=", "), "."
                           )
                     )
    if (length(xtralevs  <- setdiff(levels(object@arcs@bookkeeping[['start']]),
                                    node.labels(object)
                                    )
               )
        )
        errors  <- c(errors,
                     paste("Bookkeeping arc start nodes not listed in nodes table, e.g.",
                           paste(head(xtralevs,2), collapse=", "), "."
                           )
                     )
    if (length(xtralevs  <- setdiff(levels(object@arcs@bookkeeping[['end']]),
                                    node.labels(object)
                                    )
               )
        )
        errors  <- c(errors,
                     paste("Bookkeeping arcs' end nodes not listed in nodes table, e.g.",
                           paste(head(xtralevs,2), collapse=", "), "."
                           )
                     )
    ## Groups listed in nodes table must be same as in subproblems
    subprobs  <- unique(object@subproblems[['groups']])
    if (length(subprobs) && (length(subprobs)>1 || subprobs!=character(1)))
        {
            if (!all(unique(object@nodes[['groups']]) %in% subprobs ))
                errors  <- c(errors,
                             "Detected subproblems ('groups') in @nodes that aren't in @subproblems.")
            if (!all(unique(object@arcs@matches[['groups']]) %in% subprobs ))
                errors  <- c(errors,
                             "Detected subproblems ('groups') in @arcs@matches that aren't in @subproblems.")
            if (!all(unique(object@arcs@bookkeeping[['groups']]) %in% subprobs ))
                errors  <- c(errors,
                             "Detected subproblems ('groups') in @arcs@bookkeeping that aren't in @subproblems.")
    }
    if (length(errors)==0) TRUE else errors      
})

setClass("FullmatchMCFSolutions",
         contains="MCFSolutions",
         prototype = prototype(subproblems=new('SubProbInfo'), nodes=new('NodeInfo'),
                               arcs=new('ArcInfo'))
         )
setValidity("FullmatchMCFSolutions", function(object){
    errors  <- character(0)
    if ( nrow(object@nodes) &
         !setequal(node.labels(object)[is.na(object@nodes[['upstream_not_down']])] ,
                  c("(_Sink_)", "(_End_)") )
        )
        errors  <- c(errors,
                     "Need '(_Sink_)', '(_End_)' nodes.")
    ## (Could probably do more checking here...)
    if (length(errors)==0) TRUE else errors    
})
####################################################################
##########                  Methods            #####################
####################################################################

##* rbind.data.frame w/ revised default arguments
##*
##* @param ... data frames
##* @return data frame
##* @keywords internal
rbind_data_frame  <- function(...)
    rbind.data.frame(..., make.row.names=FALSE, stringsAsFactors=FALSE)

##* Generate combine (`c()`) functions for S4-typed data frames
##* 
##* @param thetype character, name of S4 type
##* @return function to do combine (`c()`) on object of class thetype
##* @keywords internal
typed_dataframe_combiner  <- function(thetype)
{
    function(x, ...)
    {
        ans  <- rbind_data_frame(...)
        if (!missing(x))
            ans  <- rbind_data_frame(x, ans)
        new(thetype, ans)
    }
}

setMethod("c", signature(x="SubProbInfo"),
          definition=typed_dataframe_combiner(thetype="SubProbInfo")
          )
setMethod("c", signature(x="NodeInfo"),
          definition=typed_dataframe_combiner(thetype="NodeInfo")
          )
setMethod("c", signature(x="ArcInfo"),
          definition=function(x, ...) {
              objs = list(...)
              
              mlist  <- lapply(objs, function(x) slot(x, "matches"))
              ans_m  <- do.call(rbind_data_frame, mlist)
                                
              blist  <- lapply(objs, function(x) slot(x, "bookkeeping"))
              ans_b  <- do.call(rbind_data_frame, blist)

              if (!missing(x))
              {
                  ans_m  <- rbind_data_frame(slot(x, "matches"), ans_m)
                  ans_b  <- rbind_data_frame(slot(x, "bookkeeping"), ans_b)
              }
              
              new("ArcInfo", matches=ans_m, bookkeeping=ans_b)
          })
setMethod("c", signature(x="MCFSolutions"),
          definition=function(x, ...) {
              objs  <-  list(...)
              if (!missing(x)) objs  <- c(list(x), objs)
              ans  <- new("MCFSolutions")
              ## combine NodeInfo slots first, creating a new version
              ## of the row names that's free of duplicates. Among
              ## other things, this makes the bookkeeping
              ## node levels unique by appending subproblem string. Then
              ## update factor levels in remaining slots, before
              ## attempting to combine them. 
              theslots  <- names(getSlots("MCFSolutions"))
              combined_slotvalues  <-
                  sapply(theslots, 
                         function(theslot){
                             as_list  <- lapply(objs, function(x) slot(x, theslot))
                             do.call(c, as_list)
                         }, simplify=FALSE, USE.NAMES=TRUE)
              for (theslot in theslots)
                  slot(ans, theslot)  <- combined_slotvalues[[theslot]]
              ans
          })

setMethod("c", signature(x="FullmatchMCFSolutions"),
          definition=function(x, ...) {
              ans  <- callNextMethod()
              ans  <- as(ans, "FullmatchMCFSolutions")
              ans
          } )

### node information getter:
setGeneric("nodeinfo", function(x) standardGeneric("nodeinfo"))
### (not sure yet that setter will be needed)
###setGeneric("nodeinfo<-", function(x, value) standardGeneric("nodeinfo<-"))
setMethod("nodeinfo", "NodeInfo", function(x) x)
setMethod("nodeinfo", "MCFSolutions", function(x) x@nodes)
setOldClass(c("optmatch", "factor")) # redundant given similar declaration
                                     # elsewhere; removes spurious warning.
setMethod("nodeinfo", "optmatch", function(x) {
    mcfs  <- attr(x, "MCFSolutions")
    if (is.null(mcfs)) NULL else nodeinfo(mcfs)
})
setMethod("nodeinfo", "ANY", function(x) NULL)

## node labels
setGeneric("node.labels", function(x) standardGeneric("node.labels"))
 setGeneric("node.labels<-", function(x, value) standardGeneric("node.labels<-"))
setMethod("node.labels", "NodeInfo",
          function(x) setNames(row.names(x), nm=x[['name']])
          )
setMethod("node.labels", "MCFSolutions", function(x) node.labels(x@nodes))
setMethod("node.labels", "optmatch", function(x) {
    mcfs  <- attr(x, "MCFSolutions")
    if (is.null(mcfs)) NULL else node.labels(mcfs)
})
setMethod("node.labels", "ANY", function(x) NULL)
setMethod("node.labels<-", "NodeInfo",
          function(x, value) {
              row.names(x) <- value
              x
          } 
          )
setMethod("node.labels<-", "MCFSolutions",
          function(x, value) {
              oldlabels  <- node.labels(x)
              stopifnot(length(oldlabels)==length(value))
              names(value)  <- oldlabels
              node.labels(x@nodes) <- value
              ## Now have to touch up all the factor levels!
              x@arcs@matches[['upstream']]  <-
                  factor(x@arcs@matches[['upstream']], levels=oldlabels,
                         labels = value[levels(x@arcs@matches[['upstream']])])
              x@arcs@matches[['downstream']]  <-
                  factor(x@arcs@matches[['downstream']], levels=oldlabels,
                         labels = value[levels(x@arcs@matches[['downstream']])])
              x@arcs@bookkeeping[['start']]  <-
                  factor(x@arcs@bookkeeping[['start']], levels=oldlabels,
                         labels = value[levels(x@arcs@bookkeeping[['start']])])
              x@arcs@bookkeeping[['end']]  <-
                  factor(x@arcs@bookkeeping[['end']], levels=oldlabels,
                         labels = value[levels(x@arcs@bookkeeping[['end']])]) 
              x
          } 
          )
setMethod("node.labels<-", "optmatch", function(x, value) {
    mcfs  <- attr(x, "MCFSolutions")
    if (!is.null(mcfs))
    {
        node.labels(mcfs)  <- value
        attr(x, "MCFSolutions")  <- mcfs
    }
    x
})
setMethod("node.labels<-", "ANY", function(x, value) NULL)

##'
##' This implementation does *not* drop levels of factor
##' variables (such as `groups`), in order to facilitate interpretation
##' and re-stitching together.
##' @title Filter min-cost flow information table(s) by subproblem ID
##' @param x (currently only implemented for NodeInfo objects)
##' @param ... additional parameters, including `groups`
##' @return NodeInfo
##' @author Ben Hansen
##' @keywords internal
setGeneric("filter_by_subproblem", function(x, ...) standardGeneric("filter_by_subproblem"))
setMethod("filter_by_subproblem", "NodeInfo", function(x, groups) {
    ans  <-  if (length(groups)<=1) {
                 x[x[['groups']]==groups, , drop=FALSE]
                 } else x[x[['groups']] %in% groups, , drop=FALSE]
    new("NodeInfo", ans)
    })
setMethod("filter_by_subproblem", "ANY", function(x, groups) stop("currently just implemented for NodeInfo objects"))
