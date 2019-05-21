####################################################################
########  Classes for storing information about solutions of #######
########  Min-Cost-Flow representations of matching problems #######
#######  (See vignette "MCFSolutions" for class descriptions.) #####
####################################################################
setClass("SubProbInfo", contains="data.frame",
         prototype=
             prototype(data.frame(subproblem=character(0), hashed_dist=character(0),
                              resolution=double(0), exceedance=double(0),
                              CS_orig_dist=logical(0), stringsAsFactors=FALSE)
                  )
         )
setValidity("SubProbInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:6]==
             c("subproblem","flipped", "hashed_dist","resolution","exceedance","CS_orig_dist")))
        errors  <- c(errors,
                     'Cols 1-6 should be:\n\t c("subproblem","flipped", "hashed_dist","resolution","exceedance","CS_orig_dist")')
    if (!all(vapply(object[c(1,3)], is.character, logical(1))))
        errors  <- c(errors,
                     'Cols 1,3 should have type character.')
    if (!all(vapply(object[4:5], is.double, logical(1))))
        errors  <- c(errors,
                     'Cols 4,5 should have type double.')
    if (!all(vapply(object[c(2,6)], is.logical, logical(1))))
        errors  <- c(errors,
                     'Cols 2,6 should have type logical.')
    if (anyDuplicated(object[['subproblem']]))
        errors  <- c(errors,
                     'Duplicates in "subproblem", or subproblems with same name.')
    if (length(errors)==0) TRUE else errors  
})
setClass("NodeInfo", contains="data.frame",
         prototype=
             prototype(data.frame(name=character(0), price=double(0),
                                  upstream_not_down=logical(0), supply=integer(0),
                                  subproblem=character(0), stringsAsFactors=FALSE)
                       )
         )
setValidity("NodeInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:5]==
             c("name", "price", "upstream_not_down", "supply", "subproblem")))
        errors  <- c(errors,
                     'Cols 1-5 should be:\n\t c("name", "price", "upstream_not_down", "supply", "subproblem")')
    if (!all(vapply(object[c(1,5)], is.character, logical(1))==TRUE))
        errors  <- c(errors,
                     'Cols 1,5 should have type character.')
    if (!is.double(object[['price']]))
        errors  <- c(errors,
                     'Col "price" should have type double.')
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
             prototype(matches=data.frame(subproblem=character(0), treatment=character(0), 
                                          control=character(0), stringsAsFactors=FALSE),
                       bookkeeping=data.frame(subproblem=character(0), treatment=character(0), 
                                              control=character(0), flow=integer(0),
                                              stringsAsFactors=FALSE)
                       )
         )
setValidity("ArcInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object@matches)[1:3]==
             c("subproblem", "upstream",  "downstream")))
        errors  <- c(errors,
                     '@matches cols 1-3 should be:\n\t c("subproblem", "upstream",  "downstream")')
    if (!all(vapply(object@matches, is.character, logical(1))==TRUE))
        errors  <- c(errors,
                     'All columns of @matches should have type character.')
    if (!all(colnames(object@bookkeeping)[1:4]==
             c("subproblem", "start",  "end",  "flow")))
        errors  <- c(errors,
                     '@bookkeeping cols 1-4 should be:\n\t c("subproblem", "start",  "end",  "flow")')
    if (!all(vapply(object@bookkeeping[1:3], is.character, logical(1))==TRUE))
        errors  <- c(errors,
                     '@bookkeeping cols 1-3 should have type character.')
    if (!is.integer(object@bookkeeping[['flow']]))
        errors  <- c(errors,
                     '@bookkeeping col "flow" should have type integer.')
    if (length(errors)==0) TRUE else errors      
})

setClass("MatchablesInfo", contains="data.frame",
         prototype=
             prototype(data.frame(name=character(0), 
                                  row_unit=character(0), 
                                  subproblem=character(0), stringsAsFactors=FALSE)
                       )
         )
setValidity("MatchablesInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)==
             c("name", "row_unit", "subproblem") ) )
        errors  <- c(errors,
                     'Columns should be:\n\t c("name", "row_unit", "subproblem")')
    if (!all(vapply(object[c(1,3)], is.character, logical(1))==TRUE))
        errors  <- c(errors,
                     'Cols 1,3 should have type character.')
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
    subprobs  <- unique(object@subproblems[['subproblem']])
    if (!all(unique(object@nodes[['subproblem']]) %in% subprobs ))
        errors  <- c(errors,
                     "Detected subproblems in @nodes that aren't in @subproblems.")
    if (!all(unique(object@arcs@matches[['subproblem']]) %in% subprobs ))
        errors  <- c(errors,
                     "Detected subproblems in @arcs@matches that aren't in @subproblems.")
    if (!all(unique(object@arcs@bookkeeping[['subproblem']]) %in% subprobs ))
        errors  <- c(errors,
                     "Detected subproblems in @arcs@bookkeeping that aren't in @subproblems.")
    if (!all(unique(object@matchables[['subproblem']]) %in% subprobs ))
        errors  <- c(errors,
                     "Detected subproblems in @matchables that aren't in @subproblems.")
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
