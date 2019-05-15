### See vignette "dual-MCF-solution" for descriptions of these classes.
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

setClass("ArcInfo", representation(matches="data.frame", bookkeeping="data.frame"))
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

setClass("MCFSolutions", representation(subproblems='SubProbInfo',nodes='NodeInfo',
                                        arcs='ArcInfo',matchables="MatchablesInfo"))
