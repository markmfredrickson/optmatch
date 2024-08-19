setOldClass("data.frame")
####################################################################
########  Classes for storing information about solutions of #######
########  Min-Cost-Flow representations of matching problems #######
#######  (See vignette "MCFSolutions" for class descriptions.) #####
####################################################################
setClass("SubProbInfo", contains="data.frame",
         prototype=
             prototype(data.frame(groups=character(1), flipped=NA, hashed_dist=NA_character_,
                              resolution=1, primal_value=NA_real_,
                              dual_value=NA_real_, feasible=NA, exceedance=NA_real_, stringsAsFactors=FALSE)
                  )
         )
setValidity("SubProbInfo", function(object){
    errors <- character(0)
    if (!all(colnames(object)[1:8]==
             c("groups","flipped", "hashed_dist","resolution","primal_value","dual_value", "feasible", "exceedance")))
        errors  <- c(errors,
                     'Cols 1-8 should be:\n\t c("groups","flipped", "hashed_dist","resolution","primal_value","dual_value", "feasible", "exceedance")')
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
    if (!all(colnames(object)[1:4]==
             c("name", "price", "upstream_not_down", "supply")))
        errors  <- c(errors,
                     'Cols 1-4 should be:\n\t c("name", "price", "upstream_not_down", "supply")')
    if (!is.character(object[['name']]))
        errors  <- c(errors,
                     'Col "name" should have type character.')
    if (any(colnames(object)=="groups") && !is.factor(object[['groups']]))
        errors  <- c(errors,
                     'Col "groups" should have type factor.')
    if (nlevels(factor(object[['groups']]))<=1 & # skipping this check for multi-group NodeInfo's in order
        anyDuplicated(object[['name']])          # to avoid slowing down c() methods, which (as of this
        )                                        # writing) invoke validObject() indirectly via new().
        errors <- c(errors,
                    'W/in levels of `groups`, `name`s must be unique.')
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

##' @importFrom tibble as_tibble
##' @importFrom tibble enframe
setAs("NodeInfo", "tbl_df", function(from){
    nl  <- node.labels(from)
    nl  <- factor(nl, levels=nl)
    nl  <- enframe(nl, name= NULL, value = "nodelabels")
    from  <- asS3(from) # circumvent as.tibble warning about dropping S4
    ans <- as_tibble(from, rownames=NULL)
    cbind(ans, nl)
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
    if (!identical(levels(object@matches[['upstream']]), levels(object@matches[['downstream']])))
        errors  <- c(errors, "node columns of @matches should have same levels()")
    if (!all(colnames(object@bookkeeping)[1:5]==
             c("groups", "start",  "end",  "flow", "capacity")))
        errors  <- c(errors,
                     '@bookkeeping cols 1-5 should be:\n\t c("groups", "start",  "end",  "flow", "capacity")')
    if (!all(vapply(object@bookkeeping[1:3], is.factor, logical(1))))
        errors  <- c(errors,
                     '@bookkeeping cols 1-3 should have type factor.')
    if (!identical(levels(object@bookkeeping[['start']]), levels(object@bookkeeping[['end']])))
        errors  <- c(errors,
                     "node columns of @bookkeeping should have same levels()")
    if (!identical(levels(object@matches[['upstream']]), levels(object@bookkeeping[['start']])))
        errors  <- c(errors,
                     "node columns of @matches, @bookkeeping should have same levels()")
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
    ## levels set, namely node.labels(object) (i.e. row.names(object@nodes) ).  Former
    ## requirement is enforced within the ArcInfo validity checker; here we confirm
    ## that this level set matches what's in the nodes table.
    if (!identical(row.names(nodeinfo(object)),
                   levels(object@arcs@matches[['upstream']])
                   )
        )
        {
            if (length(xtralevs  <- setdiff(node.labels(object),
                                            levels(object@arcs@matches[['upstream']])
                                            )
                       )
                )
                errors  <- c(errors,
                             paste("@nodes lists nodes not in the levels of @arcs's nodes columns, e.g.",
                                   paste(head(xtralevs,2), collapse=", "), "."
                                   )
                             )
            if (length(xtralevs  <- setdiff(levels(object@arcs@matches[['upstream']]),
                                            node.labels(object)
                                            )
                       )
                )
                errors  <- c(errors,
                             paste("@arcs's nodes columns have levels not appearing in @nodes",
                                   paste(head(xtralevs,2), collapse=", "), "."
                                   )
                             )
        }
    ## Nodes table must have groups column
    if (!any(colnames(object@nodes)=="groups"))
        errors  <- c(errors,
                     "nodes table must have 'groups' column")
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

##' Combine objects
##' @param x object of particular class
##' @param ... Various objects
##' @return Combined objects
##' @rdname mcf_c_fns
##' @export
setMethod("c", signature(x="SubProbInfo"),
          definition=typed_dataframe_combiner(thetype="SubProbInfo")
          )
##' @rdname mcf_c_fns
##' @export
setMethod("c", signature(x="NodeInfo"),
          definition=typed_dataframe_combiner(thetype="NodeInfo")
          )
##' @rdname mcf_c_fns
##' @export
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
##' @rdname mcf_c_fns
##' @export
setMethod("c", signature(x="MCFSolutions"),
          definition=function(x, ...) {
              objs  <-  list(...)
              if (!missing(x)) objs  <- c(list(x), objs)

              ans  <- new("MCFSolutions")
              methods::slot(ans, "subproblems")  <- do.call(c, lapply(objs, slot, "subproblems"))
              methods::slot(ans, "nodes")  <- do.call(c, lapply(objs, slot, "nodes"))

              ## This creates new node.labels (that are free of duplicates).
              new_nl  <- node.labels(slot(ans, "nodes"))
              ## Propagating these new labels into the ArcInfo's takes a few steps.

              n_nodes_by_obj  <- vapply(objs, function(obj) nrow(nodeinfo(obj)), integer(1))
              offsets  <- if (length(objs)>1) c(0L, cumsum(n_nodes_by_obj)[-length(objs)]) else 0L
              for (ii in 1L: length(objs))
                 methods::slot(objs[[ii]], "arcs")  <-
                      revise_ArcInfo_nodelabels(slot(objs[[ii]], "arcs"),
                                                new_nl,
                                                offsets[ii]+1L:n_nodes_by_obj[ii])

              methods::slot(ans, "arcs")  <- do.call(c, lapply(objs, slot, "arcs"))
              ans
          })

##' @rdname mcf_c_fns
##' @export
setMethod("c", signature(x="FullmatchMCFSolutions"),
          definition=function(x, ...) {
              ans  <- methods::callNextMethod()
              ans  <- as(ans, "FullmatchMCFSolutions")
              ans
          } )

### node information getter:
setGeneric("nodeinfo", function(x) standardGeneric("nodeinfo"))
setMethod("nodeinfo", "NodeInfo", function(x) x)
setMethod("nodeinfo", "MCFSolutions", function(x) x@nodes)
setOldClass(c("optmatch", "factor")) # redundant given similar declaration
                                     # elsewhere; removes spurious warning.
setMethod("nodeinfo", "optmatch", function(x) {
    mcfs  <- attr(x, "MCFSolutions")
    if (is.null(mcfs)) NULL else nodeinfo(mcfs)
})
setMethod("nodeinfo", "ANY", function(x) NULL)
## node information setter.  Note presumptions of a single `group`,
## and that the new nodes are a superset of the old ones. (The one-
## group presumption is only b/c I didn't need it for more than 1 group.)
## Note also that they're lined up on the basis of their name column,
## not their node.labels.
setGeneric("nodeinfo<-", function(x, value) standardGeneric("nodeinfo<-"))
setMethod("nodeinfo<-", c(x="MCFSolutions", value="NodeInfo"),
          function(x, value) {
              stopifnot(isTRUE(methods::validObject(x)),
                        length(unique(nodeinfo(x)[['groups']]))<=1,
                        !any(is.na(positions  <-
                                       match(nodeinfo(x)[['name']], value[['name']])
                                   )
                             )
                        )
              ## prepare to update the various factor levels
              oldlevs  <- node.labels(x)
              newlevs  <- node.labels(value)
              intermediatelevs  <- newlevs
              intermediatelevs[positions]  <- oldlevs
              ## first the easy part:
              x@nodes  <- value
              ## now we cycle through and update all the levels
              x@arcs@matches[['upstream']]  <- factor(x@arcs@matches[['upstream']],
                                                 levels=intermediatelevs,
                                                 labels=newlevs)
              x@arcs@matches[['downstream']]  <- factor(x@arcs@matches[['downstream']],
                                                 levels=intermediatelevs,
                                                 labels=newlevs)
              x@arcs@bookkeeping[['start']]  <- factor(x@arcs@bookkeeping[['start']],
                                                 levels=intermediatelevs,
                                                 labels=newlevs)
              x@arcs@bookkeeping[['end']]  <- factor(x@arcs@bookkeeping[['end']],
                                                 levels=intermediatelevs,
                                                 labels=newlevs)
              x
})
setMethod("nodeinfo<-", "ANY", function(x, value) stop("Not implemented."))

## node labels
setGeneric("node.labels", function(x) standardGeneric("node.labels"))
setGeneric("node.labels<-", function(x, value) standardGeneric("node.labels<-"))
#' @importFrom stats setNames
setMethod("node.labels", "NodeInfo",
          function(x) stats::setNames(row.names(x), nm=x[['name']])
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
##' @title Reset implicit node labels of an ArcInfo object
##' @param x an ArcInfo object
##' @param new character; the new node labels (level sets for factors encoding arc start or end nodes)
##' @param old_positions integer; positions for old levels with new levels vector
##' @return ArcInfo object with new levels
##' @author Ben Hansen
##' @keywords internal
revise_ArcInfo_nodelabels  <- function(x, new, old_positions=1L:length(new))
{
    stopifnot(is(x, "ArcInfo"),
              validObject(x), # maybe remove after development?
              is.character(new),
              is.integer(old_positions),
              all(1L<=old_positions),
              all(old_positions<=length(new)),
              !any(duplicated(old_positions)),
              length(old_positions)==nlevels(x@matches[['upstream']])
              )
    n_old  <- length(old_positions)
    n_new  <- length(new)
    oldlevs  <- levels(x@matches[['upstream']])
    ## we need a temporary level set that extends the old level set
    ## as necessary in order to match the length of `new`.
    oldlevs_padded  <- character(n_new)
    oldlevs_padded[old_positions]  <- oldlevs
    ## If we're not only renaming but also adding levels, we'll
    ## need to pad `oldlevs`.  It doesn't matter what the padding is
    ## so long as it's distinct from the actual old levels.
    if (n_old < n_new)
        oldlevs_padded[-old_positions]  <-
            make.unique(c(oldlevs, new[-old_positions])
                        )[(n_old+1L):n_new]
    ## Finally, touch up all the factor levels.
    x@matches[['upstream']]  <-
        factor(x@matches[['upstream']],
               levels=oldlevs_padded,
               labels = new)
    x@matches[['downstream']]  <-
        factor(x@matches[['downstream']],
               levels=oldlevs_padded,
               labels = new)
    x@bookkeeping[['start']]  <-
        factor(x@bookkeeping[['start']],
               levels=oldlevs_padded,
               labels = new)
    x@bookkeeping[['end']]  <-
        factor(x@bookkeeping[['end']],
               levels=oldlevs_padded,
               labels = new)
    x
    }
setMethod("node.labels<-", "MCFSolutions",
          function(x, value) {
              oldlabels  <- node.labels(x)
              stopifnot(length(oldlabels)==length(value))
              names(value)  <- oldlabels
              node.labels(x@nodes) <- value
              x@arcs  <-
                  revise_ArcInfo_nodelabels(x@arcs, value)
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

## Function to update a node table's prices and supplies
## using from a second node table carrying price and supply
## info about a subset of the first table's nodes. The two
## tables will be aligned using their "name" columns,
## not their node.labels (rownames).
update.NodeInfo  <- function(object, new, ...)
{
    stopifnot(is(object, "NodeInfo"), is(new, "NodeInfo"),
              !any(is.na(positions  <- match(new[['name']], object[['name']])))
              )
    price_col_position_old  <- which(object@names=="price")
    price_col_position_new  <- which(new@names=="price")
    object@.Data[[price_col_position_old]][positions]  <-
        new@.Data[[price_col_position_new]]

    supply_col_position_old  <- which(object@names=="supply")
    supply_col_position_new  <- which(new@names=="supply")
    object@.Data[[supply_col_position_old]][positions]  <-
        new@.Data[[supply_col_position_new]]

    object
}
##* `dplyr::filter()` method for NodeInfo's.
filter.NodeInfo  <- function(.data, ...) {
    x  <- as(.data, "tbl_df")
    ans  <- filter(x, ...)
    nl  <- which(colnames(ans)=="nodelabels")
    ans <- as.data.frame(ans[-nl],
                         row.names=as.character(ans[[nl]])
                         )
    new("NodeInfo", ans)
}
##*
##* This implementation does *not* drop levels of factor
##* variables (such as `groups`), in order to facilitate interpretation
##* and re-stitching together.
##* @title Filter min-cost flow information table(s) by subproblem ID
##* @param x (currently only implemented for NodeInfo objects)
##* @param ... additional parameters, including `groups`
##* @return NodeInfo
##* @author Ben Hansen
##* @keywords internal
setGeneric("filter_by_subproblem", function(x, ...) standardGeneric("filter_by_subproblem"))
setMethod("filter_by_subproblem", "NodeInfo", function(x, groups) {
    ans  <-  if (length(groups)<=1) {
                 x[x[['groups']]==groups, , drop=FALSE]
                 } else x[x[['groups']] %in% groups, , drop=FALSE]
    new("NodeInfo", ans)
    })
setMethod("filter_by_subproblem", "ANY", function(x, groups) stop("currently just implemented for NodeInfo objects"))
