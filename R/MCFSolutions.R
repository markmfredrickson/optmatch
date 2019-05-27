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
                              dual_value=NA_real_, feasible=NA, stringsAsFactors=FALSE)
                  )
         )
setValidity("SubProbInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:7]==
             c("groups","flipped", "hashed_dist","resolution","lagrangian_value","dual_value", "feasible")))
        errors  <- c(errors,
                     'Cols 1-7 should be:\n\t c("groups","flipped", "hashed_dist","resolution","lagrangian_value","dual_value", "feasible")')
    if (!all(vapply(object[c(1,3)], is.character, logical(1))))
        errors  <- c(errors,
                     'Cols 1,3 should have type character.')
    if (!all(vapply(object[4:6], is.double, logical(1))))
        errors  <- c(errors,
                     'Cols 4-6 should have type double.')
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
                     'in @bookkeeping, flow can be now greater than capacity.')
    if (length(errors)==0) TRUE else errors      
})

setClass("MatchablesInfo", contains="data.frame",
         prototype=
             prototype(data.frame(name = character(0),
                                  row_unit = logical(0),
                                  groups = factor(), stringsAsFactors=FALSE)
                       )
         )
setValidity("MatchablesInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)==
             c("name", "row_unit", "groups") ) )
        errors  <- c(errors,
                     'Columns should be:\n\t c("name", "row_unit", "groups")')
    if (!is.character(object[['name']]))
        errors  <- c(errors,
                     'Cols "name" should have type character.')
    if (!is.factor(object[['groups']]))
        errors  <- c(errors,
                     'Cols "groups" should have type factor.')
    if ( !is.logical(object[['row_unit']]) )
        errors  <- c(errors,
                     "'row_unit' col should have type logical.")
    
    if (length(errors)==0) TRUE else errors  
})

setClass("MCFSolutions", slots=c(subproblems='SubProbInfo',nodes='NodeInfo',
                                 arcs='ArcInfo',matchables="MatchablesInfo"),
         prototype = prototype(subproblems=new('SubProbInfo'), nodes=new('NodeInfo'),
                               arcs=new('ArcInfo'),matchables=new("MatchablesInfo"))
         )
setValidity("MCFSolutions", function(object){
    errors  <- character(0)
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
            if (!all(unique(object@matchables[['groups']]) %in% subprobs ))
                errors  <- c(errors,
                             "Detected subproblems ('groups') in @matchables that aren't in @subproblems.")
    }
    if (length(errors)==0) TRUE else errors      
})

setClass("FullmatchMCFSolutions",
         contains="MCFSolutions",
         prototype = prototype(subproblems=new('SubProbInfo'), nodes=new('NodeInfo'),
                               arcs=new('ArcInfo'),matchables=new("MatchablesInfo"))
         )
setValidity("FullmatchMCFSolutions", function(object){
    errors  <- character(0)
    if ( nrow(object@nodes) &
         !setequal(object@nodes[['name']][is.na(object@nodes[['upstream_not_down']])] ,
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
setMethod("c", signature(x="MatchablesInfo"),
          definition=typed_dataframe_combiner(thetype="MatchablesInfo")
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
