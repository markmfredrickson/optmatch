####################################################################
########  Classes for storing information about solutions of #######
########  Min-Cost-Flow representations of matching problems #######
#######  (See vignette "MCFSolutions" for class descriptions.) #####
####################################################################
setClass("SubProbInfo", contains="data.frame")
setValidity("SubProbInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:5]==
             c("subproblem","hashed_dist","resolution","exceedance","CS_orig_dist")))
        errors  <- c(errors,
                     'Cols 1-5 should be:\n\t c("subproblem","hashed_dist","resolution","exceedance","CS_orig_dist")')
    if (!all(sapply(object[1:2], is.character)==TRUE))
        errors  <- c(errors,
                     'Cols 1,2 should have type character.')
    if (!all(sapply(object[3:4], is.numeric)==TRUE))
        errors  <- c(errors,
                     'Cols 3,4 should have type numeric.')
    if (!is.logical(object[[5]]))
        errors  <- c(errors,
                     'Col 5 should have type logical.')
    if (length(errors)==0) TRUE else errors  
})
setClass("NodeInfo", contains="data.frame")
setValidity("NodeInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:5]==
             c("name", "price", "kind", "supply", "subproblem")))
        errors  <- c(errors,
                     'Cols 1-5 should be:\n\t c("name", "price", "kind", "supply", "subproblem")')
    if (!all(sapply(object[c(1,3,5)], is.character)==TRUE))
        errors  <- c(errors,
                     'Cols 1,3,5 should have type character.')
    if (!is.double(object[['price']]))
        errors  <- c(errors,
                     'Col "price" should have type double.')
    if (!is.integer(object[['supply']]))
        errors  <- c(errors,
                     'Col "supply" should have type integer.')    
    if ( !all(object[['kind']] %in% c("treatment", "control", "bookkeeping")) )
        errors  <- c(errors,
                     "'kind' values other than 'treatment', 'control' or 'bookkeeping'.")
    if (length(errors)==0) TRUE else errors      
})

setClass("ArcInfo", slots=c(matches="data.frame", bookkeeping="data.frame"))
setValidity("ArcInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object@matches)[1:3]==
             c("subproblem", "treatment",  "control")))
        errors  <- c(errors,
                     '@matches cols 1-3 should be:\n\t c("subproblem", "treatment",  "control")')
### <ADD COLUMN TYPE CHECKS HERE>
    
    if (!all(colnames(object@bookkeeping)[1:4]==
             c("subproblem", "startnode",  "endnode",  "flow")))
        errors  <- c(errors,
                     '@bookkeeping cols 1-4 should be:\n\t c("subproblem", "startnode",  "endnode",  "flow")')
    if (!is.integer(object@bookkeeping[['flow']]))
        errors  <- c(errors,
                     '@bookkeeping col "flow" should have type integer.')
### <ADD FURTHER COLUMN TYPE CHECKS HERE>
    if (length(errors)==0) TRUE else errors      
})

setClass("MatchablesInfo", contains="data.frame")
setValidity("MatchablesInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)==
             c("name", "kind", "subproblem") ) )
        errors  <- c(errors,
                     'Columns should be:\n\t c("name", "kind", "subproblem")')
    if (!all(sapply(object, is.character)==TRUE))
        errors  <- c(errors,
                     'All columns should have type character.')
    if ( !all(object[['kind']] %in% c("treatment", "control")) )
        errors  <- c(errors,
                     "'kind' values other than 'treatment' or 'control'.")
    
    if (length(errors)==0) TRUE else errors  
})

setClass("MCFSolutions", slots=c(subproblems='SubProbInfo',nodes='NodeInfo',
                                        arcs='ArcInfo',matchables="MatchablesInfo"))
setValidity("MCFSolutions", function(object){
    errors  <- character(0)
    if (any(duplicated(object@subproblems[['subproblem']])))
        errors  <- c(errors,
                     "@subproblems lists duplicates or subproblems with same name.")    
    if (length(errors)==0) TRUE else errors      
})

##' Check cross-compatibility of slots constituting
##' an MCFSolutions object. 
##'
##' These checks are potentially time-consuming, so they're
##' parked here instead of in officialy validity-checker
##' in order to avoid replicating them when (previously checked)
##' MDFSolutions objects are being combined with `c()`. 
##' @title Additional validity check for MCFSolutions
##' @param object MCFSolutions
##' @return TRUE if valid, or character vector of errors otherwise
##' @keywords internal
validMCFSolutions  <- function(object) {
    stopifnot("MCFSolutions" %in% is(object))
    errors  <- validObject(object)
    if (isTRUE(errors)) errors  <- character(0)
    
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

    if (length(errors)==0) TRUE else paste(errors, "\n")
    }
####################################################################
##########                  Methods            #####################
####################################################################

#' rbind.data.frame w/ revised default arguments
#'
#' @param ... data frames
#' @return data frame
#' @keywords internal
rbind_data_frame  <- function(...)
    rbind.data.frame(..., make.row.names=FALSE, stringsAsFactors=FALSE)

#' Generate combine (`c()`) functions for S4-typed data frames
#' 
#' @param thetype character, name of S4 type
#' @return function to do combine (`c()`) on object of class thetype
#' @keywords internal
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
              theslots  <- names(getSlots("MCFSolutions"))
              combined_slotvalues  <-
                  sapply(theslots, 
                         function(theslot){
                             as_list  <- lapply(objs, function(x) slot(x, theslot))
                             do.call(c, as_list)
                         }, USE.NAMES=TRUE)
              do.call(new, c(list(Class="MCFSolutions"), combined_slotvalues))
          })
