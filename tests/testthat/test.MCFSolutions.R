context("MCFSolutions & co: S4 classes to encode min-cost-flow solutions")

test_that("Instantiation & validity", {
    expect_true(validObject(new("SubProbInfo"))) # test the prototype
    expect_silent(spi1  <- new("SubProbInfo", data.frame(groups=c('a','b'), flipped=logical(2), hashed_dist=c('a','b'),
                                           resolution=c(1,10), lagrangian_value=c(.5, 2), dual_value=c(0,1),
                                           feasible=c(TRUE, FALSE), exceedance=1, stringsAsFactors=F)
                               )
                  )
    expect_true(validObject(spi1))

    spi2  <- spi1
    colnames(spi2)[1]  <- "Subprob"
    expect_error(validObject(spi2), "Cols 1-8 should be")

    expect_true(validObject(new("NodeInfo"))) # test the prototype    
    expect_silent(ni1  <- new("NodeInfo",
                             data.frame(name='a', price=0.5, upstream_not_down=TRUE,
                                        supply=1L, groups = as.factor('b'),
                                        stringsAsFactors=F)
                             )
                  )
    expect_silent(ni1f  <- new("NodeInfo",
                               data.frame(name=c('b', 'c', 'd',
                                                 '(_Sink_)', '(_End_)'),
                                          price=0.5,
                                          upstream_not_down=c(TRUE, FALSE,
                                                              FALSE, NA, NA),
                                         supply=c(1L,0L,0L,-1L,-2L), groups = as.factor('b'),
                                         stringsAsFactors=F)
                              )
                  )
    expect_equivalent(node.labels(ni1f), as.character(1:5))
    expect_named(node.labels(ni1f),
                 c('b', 'c', 'd', '(_Sink_)', '(_End_)')
                 )
    expect_silent(node.labels(ni1f)  <- ni1f[['name']])
    expect_equivalent(node.labels(ni1f), ni1f[['name']])
    expect_error(new("NodeInfo",
                      data.frame(name='a', price='foo', upstream_not_down=FALSE,
                                 supply=1L, groups = as.factor('b'),
                                 stringsAsFactors=F)
                     ),
                 "should be a numeric" # Not sure it's necessary, but insisting
                  )                        # that 'price' be double not integer
    expect_error(new("NodeInfo",
                      data.frame(name=rep('a', 2), price=0, upstream_not_down=FALSE,
                                 supply=1L, groups = as.factor('b'),
                                 stringsAsFactors=F)
                     ),
                 "unique" 
                  )
    expect_true(validObject(new("ArcInfo"))) # test the prototype
    expect_silent(ai  <- new("ArcInfo",
                             matches=data.frame(groups = as.factor('a'), upstream = factor('b', levels=node.labels(ni1f)),
                                                downstream = factor(c('c','d'), levels=node.labels(ni1f)),
                                                stringsAsFactors=F),
                             bookkeeping=data.frame(groups = as.factor('a'),
                                                    start = factor(c('c','d'),
                                                                   levels=node.labels(ni1f)
                                                                   ),
                                                    end = factor('(_Sink_)',
                                                                 levels=node.labels(ni1f)
                                                                 ),
                                                    flow=1L, capacity=1L,
                                                    stringsAsFactors=F)
                             )
                  )
    expect_error(new("ArcInfo",
                     matches=data.frame(groups = as.factor('a'), upstream = as.factor('b'),
                                        downstream = as.factor(c('c','d')),stringsAsFactors=F),
                     bookkeeping=data.frame(groups = as.factor('a'), start = as.factor(c('c','d')),
                                            end = as.factor('(_Sink_)'), flow=1.5,
                                            capacity=1L, stringsAsFactors=F)
                     ), "should have type integer" # Not sure it's necessary, but insisting 
                  )                                # that 'flow' be integer not double
    expect_error(new("ArcInfo",
                     matches=data.frame(groups = as.factor('a'), upstream = as.factor('b'),
                                        downstream = as.factor(c('c','d')),stringsAsFactors=F),
                     bookkeeping=data.frame(groups = as.factor('a'), start = as.factor(c('c','d')),
                                            end = as.factor('(_Sink_)'), flow=-1L,
                                            capacity=1L, stringsAsFactors=F)
                     ), "should be nonnegative" 
                  )                                
    expect_error(new("ArcInfo",
                     matches=data.frame(groups = as.factor('a'),
                                        upstream = as.factor('b'),
                                        downstream=as.factor(c('c','d')),
                                        stringsAsFactors=F),
                     bookkeeping=data.frame(groups=as.factor('a'),
                                            start=as.factor(c('c','d')),
                                            end=as.factor('(_Sink_)'), flow=2L,
                                            capacity=1L, stringsAsFactors=F)
                     ), "flow can be no greater than capacity" 
                  )                                
    
    expect_silent(mcf1f  <- new("MCFSolutions", subproblems=spi1,nodes=ni1f,arcs=ai))
    expect_silent(as(mcf1f, "FullmatchMCFSolutions"))
    expect_equal(node.labels(mcf1f), node.labels(ni1f))
    expect_silent(node.labels(mcf1f) <- paste0(node.labels(ni1f),"_") )
    expect_equivalent(node.labels(mcf1f), paste0(node.labels(ni1f), "_") )
    expect_equivalent(node.labels(mcf1f@nodes), paste0(node.labels(ni1f), "_") )
})

test_that("c() methods", {
    spi1  <- new("SubProbInfo",
                 data.frame(groups=c('a','b'), flipped=logical(2), hashed_dist=c('a','b'),
                            resolution=c(1,10), lagrangian_value=c(.5, 2), dual_value=c(0,1),
                            feasible=c(TRUE, FALSE), exceedance=1, stringsAsFactors=F)
                 )
    spi2  <- new("SubProbInfo",
                 data.frame(groups=c('c'), flipped=logical(1), hashed_dist=c('a'),
                            resolution=c(1), lagrangian_value=c(.5), dual_value=c(0),
                            feasible=c(TRUE), exceedance=1, stringsAsFactors=F)
                 )
    spi3  <- new("SubProbInfo",
                 data.frame(groups=c('d'), flipped=logical(1), hashed_dist=c('a'),
                            resolution=c(1), lagrangian_value=c(.5), dual_value=c(0),
                            feasible=c(TRUE), exceedance=1, stringsAsFactors=F)
                 )
    
    expect_silent(c(spi1, spi2))
    expect_silent(c(spi1, spi2, spi3))
    expect_silent(c(a=spi1, b=spi2)) # no confusion just b/c no `x=` arg!
    
    ni1f  <- new("NodeInfo",
                 data.frame(name=c('b', 'c', 'd',
                                   '(_Sink_)', '(_End_)'),
                            price=c(0.5, 0.5,
                                    NA_real_, # permissible for downstream nodes
                                    0.5, 0.5),
                            upstream_not_down=c(TRUE, FALSE,
                                                FALSE, NA, NA),
                            supply=c(1L,0L,0L,-1L,-2L), groups = as.factor('b'),
                            stringsAsFactors=F)
                 )
    node.labels(ni1f) <- ni1f[['name']]
    ni1f.a  <- ni1f.b <- ni1f.c  <- ni1f
    ni1f.a[,'groups']  <- factor(rep('a', nrow(ni1f)))
    ni1f.c[,'groups']  <- factor(rep('c', nrow(ni1f)))    
    expect_silent(c(ni1f.a, ni1f.b))
    expect_silent(c(ni1f.a, ni1f.b, ni1f.c))
    ni2  <- new("NodeInfo",
                data.frame(name=c(letters[2:5], '(_Sink_)', '(_End_)'), price=0.5,
                           upstream_not_down=c(TRUE, rep(FALSE,3), NA, NA),
                          supply=c(1L, rep(0L,3),-1L,-2L), groups=as.factor('c'),
                          stringsAsFactors=F)
                )
    node.labels(ni2) <- ni2[['name']]
    ni1ni2 <- c(ni1f, ni2)
    expect_equal(ni1ni2$name, c("b", "c", "d", "(_Sink_)", "(_End_)", letters[2:5], "(_Sink_)", "(_End_)"))
    expect_equal(levels(ni1ni2$groups), c("b", "c"))
    expect_named(node.labels(ni1ni2), c("b", "c", "d", "(_Sink_)", "(_End_)", letters[2:5], "(_Sink_)", "(_End_)") )
    expect_false( any(duplicated(node.labels(ni1ni2))) )

    ai1  <- new("ArcInfo",
               matches=data.frame(groups = factor('a'), upstream = factor('b'),
                                  downstream = factor(c('c','d')), stringsAsFactors=F),
               bookkeeping=data.frame(groups= factor('a'), start = factor(c('c','d')),
                                      end= factor('(_Sink_)'), flow=1L,
                                      capacity=1L, stringsAsFactors=F)
               )
    expect_silent(c(ai1, ai1))
    expect_silent(c(x=ai1, y=ai1, z=ai1))
    ai2  <- new("ArcInfo",
               matches=data.frame(groups = factor('c'), upstream = factor('b'),
                                  downstream = factor(c('c','d', 'e')), stringsAsFactors=F),
               bookkeeping=data.frame(groups = factor('c'), start = factor(c('c','d')),
                                      end = factor('(_Sink_)'), flow=1L,
                                      capacity=1L, stringsAsFactors=F)
               )

    ai1ai2 <- c(ai1, ai2)
    expect_equal(levels(ai1ai2@matches$groups), c("a", "c"))
    expect_equal(levels(ai1ai2@matches$upstream), "b")
    expect_equal(levels(ai1ai2@matches$downstream), c("c", "d", "e"))
    expect_equal(levels(ai1ai2@bookkeeping$end), "(_Sink_)")


    mcf1  <- new("MCFSolutions", subproblems=spi1, nodes=ni1f, arcs=ai1)
    expect_error(c(mcf1, mcf1), "uplicates")


    mcf2 <- new("MCFSolutions", subproblems=spi2,nodes=ni2,arcs=ai2)
    
    expect_silent(c(mcf1, mcf2))
    expect_silent(c(y=mcf1, z=mcf2))

    mcf2f  <- as(mcf2, "FullmatchMCFSolutions")
    expect_is(c(mcf2f, mcf1), "FullmatchMCFSolutions")
    expect_is(c(mcf1, mcf2f), "MCFSolutions")
})
test_that("nodeinfo getter",{
    expect_silent(mcf  <-  new("MCFSolutions")) #prelim-
    expect_true(validObject(mcf, complete=TRUE))#inaries

    expect_is(nodeinfo(mcf@nodes), "NodeInfo")
    expect_is(nodeinfo(mcf), "NodeInfo")

    data <- data.frame(z = c(rep(0,10), rep(1,5)),
                       x = rnorm(15), fac=rep(c(rep("a",2), rep("b",3)),3))
    f1 <- fullmatch(z~x, min.c=1, max.c=1, omit.fraction=.5, data = data)
    expect_is(f1, "optmatch")
    expect_false(is.null(attr(f1, "MCFSolutions")))
    expect_is(nodeinfo(f1), "NodeInfo")
    expect_null(nodeinfo(10))
})
test_that("NodeInfo to tibble converter", {
    ni1f  <- new("NodeInfo",
                 data.frame(name=c('b', 'c', 'd',
                                   '(_Sink_)', '(_End_)'),
                            price=0.5,
                            upstream_not_down=c(TRUE, FALSE,
                                                FALSE, NA, NA),
                            supply=c(1L,0L,0L,-1L,-2L), groups = as.factor('b'),
                            stringsAsFactors=F)
                 )
    expect_silent(ni_tbl  <- as(ni1f, "tbl_df"))
    expect_is(ni_tbl$nodelabels, "factor")
    expect_equivalent(as.character(ni_tbl$nodelabels),
                      as.character(1:5))
    expect_null(names(ni_tbl$nodelabels))
    node.labels(ni1f) <- ni1f[['name']]
    expect_silent(ni_tbl  <- as(ni1f, "tbl_df"))
    expect_is(ni_tbl$nodelabels, "factor")
    expect_equivalent(as.character(ni_tbl$nodelabels), 
                      c('b', 'c', 'd', '(_Sink_)', '(_End_)')
                      )
    expect_equivalent(levels(ni_tbl$nodelabels), 
                      c('b', 'c', 'd', '(_Sink_)', '(_End_)')
                      ) # default encoding would start w/ "(_End_)", "(_Sink_)"
})
test_that("Preserve levels when filtering a node info tibble",{
    ni1f  <- new("NodeInfo",
                 data.frame(name=c('b', 'c', 'd',
                                   '(_Sink_)', '(_End_)'),
                            price=0.5,
                            upstream_not_down=c(TRUE, FALSE,
                                                FALSE, NA, NA),
                            supply=c(1L,0L,0L,-1L,-2L), groups = as.factor('b'),
                            stringsAsFactors=F)
                 )
    node.labels(ni1f) <- ni1f[['name']]
    ni_tbl  <- as(ni1f, "tbl_df")
    expect_silent(ni_tbl_s  <- dplyr::filter(ni_tbl, name %in% letters))
    expect_is(ni_tbl_s$nodelabels, "factor")
    expect_equivalent(levels(ni_tbl_s$nodelabels), 
                      c('b', 'c', 'd', '(_Sink_)', '(_End_)')
                      ) 
})
test_that("Node labels getter",{
    expect_silent(mcf  <-  new("MCFSolutions")) #prelim-
    expect_true(validObject(mcf, complete=TRUE))#inaries

    expect_is(node.labels(mcf@nodes), "character")
    expect_is(node.labels(mcf), "character")

    data <- data.frame(z = c(rep(0,10), rep(1,5)),
                       x = rnorm(15), fac=rep(c(rep("a",2), rep("b",3)),3))
    f1 <- fullmatch(z~x, min.c=1, max.c=1, omit.fraction=.5, data = data)
    expect_is(f1, "optmatch")
    expect_false(is.null(attr(f1, "MCFSolutions")))
    expect_is(node.labels(f1), "character")
    expect_false(is.null(names(node.labels(f1))))
    expect_null(nodeinfo(10))

})
test_that("filtering on groups/subproblem field", {

    spi1  <- new("SubProbInfo",
                 data.frame(groups=c('a','b'), flipped=logical(2), hashed_dist=c('a','b'),
                            resolution=c(1,10), lagrangian_value=c(.5, 2), dual_value=c(0,1),
                            feasible=c(TRUE, FALSE), exceedance=1, stringsAsFactors=F)
                 )
    expect_error(filter_by_subproblem(spi1, groups="a"), "implemented")

    ni1f  <- new("NodeInfo",
                 data.frame(name=c('b', 'c', 'd',
                                   '(_Sink_)', '(_End_)'),
                            price=0.5,
                            upstream_not_down=c(TRUE, FALSE,
                                                FALSE, NA, NA),
                            supply=c(1L,0L,0L,-1L,-2L), groups = as.factor('b'),
                            stringsAsFactors=F
                            )
                 )
    node.labels(ni1f)  <- ni1f[['name']]
    expect_silent(ni1a  <- filter_by_subproblem(ni1f, groups="b"))
    expect_identical(ni1f, ni1a)
    expect_silent(ni10  <- filter_by_subproblem(ni1f, groups="a"))
    expect_is(ni10, "NodeInfo")
    expect_equal(nrow(ni10), 0L)

    ni2  <- new("NodeInfo",
                data.frame(name=c('a', '(_Sink_)', '(_End_)'), price=0.5,
                           upstream_not_down=c(FALSE, NA, NA),
                          supply=c(0L,-1L,-2L), groups=as.factor('c'),
                          stringsAsFactors=F)
               )
    ni12  <- c(ni1f, ni2)
    expect_silent(ni1b  <- filter_by_subproblem(ni12, groups="b"))
    expect_is(ni1b, "NodeInfo")
    expect_equal(nrow(ni1b), 5L)
    expect_silent(ni12a  <- filter_by_subproblem(ni12, groups=c("b","c")))
    expect_is(ni12a, "NodeInfo")
    expect_equal(nrow(ni12a), 8L)
    
    ai1  <- new("ArcInfo",
               matches=data.frame(groups = factor('a'), upstream = factor('b'),
                                  downstream = factor(c('c','d')), stringsAsFactors=F),
               bookkeeping=data.frame(groups= factor('a'), start = factor(c('c','d')),
                                      end= factor('(_Sink_)'), flow=1L,
                                      capacity=1L, stringsAsFactors=F)
               )
    expect_error(filter_by_subproblem(ai1, groups="a"), "implemented")

    mcf1  <- new("MCFSolutions", subproblems=spi1, nodes=ni1f, arcs=ai1)
    expect_error(filter_by_subproblem(mcf1, groups="a"), "implemented")
})

test_that("Potentially unusual requirements of base functions",{
    ## de-duplication of row names when stacking data.frame-s
    expect_true(formals(base::rbind.data.frame)[['make.row.names']])
    df1 <- data.frame(x=1:2, y=3:4, row.names=c('a','b'))
    df2 <- data.frame(x=3:4, y=3:4, row.names=c('a','B'))
    df3 <- data.frame(x=5:6, y=3:4, row.names=c('b','a'))
    expect_equal(row.names(rbind(df1, df2)), c("a",  "b",  "a1", "B"))
    expect_equal(row.names(rbind(df1, df3)), c("a",  "b",  "b1", "a1"))
    ## char vecs can have repeats in row names
    cvec  <- letters
    names(cvec)  <- LETTERS
    names(cvec)[1]  <- "B"
    expect_is(cvec, "character")
    expect_named(cvec, c("B", LETTERS[-1L]))
    expect_true(any(duplicated(names(cvec))))
})

test_that("NodeInfo subsetting", {
    expect_silent(ni0_o  <- nodes_shell_fmatch(c(1,2), c(3,4)))
    expect_is(ni0_o, "NodeInfo")
    expect_equal(nrow(ni0_o), 6)
    expect_silent(ni0_n  <- filter(ni0_o, name!=2 & name!=4))
    expect_is(ni0_n, "NodeInfo")
    expect_equal(nrow(ni0_n), 4)    
})
test_that("Pull updated prices & supplies into a NodeInfo",{
    expect_silent(ni0_o  <- nodes_shell_fmatch(c(1,2), c(3,4)))
    expect_true(all(ni0_o[['price']]==0))
    ni0_n  <- filter(ni0_o, name!=2 & name!=4)
    expect_equal(nrow(ni0_n), 4) # unimportant in itself but presumed by next few lines
    ni0_n@.Data[[which(ni0_n@names=="price")]]  <- rep(1, 4)
    ni0_n@.Data[[which(ni0_n@names=="supply")]]  <- c(2, 0, -1, -1)
    ## (Now ni0_n has the form of an actual price/supply combo)
    expect_silent(ni0  <- update.NodeInfo(ni0_o, ni0_n))
    expect_equal(ni0[['price']], c(1, 0, 1, 0, 1, 1))
    expect_equal(ni0[['supply']], c(2, 0, 0, 0, -1, -1))

    expect_equal(update(ni0_o, ni0_n), ni0) # confirm `update()` dispatch
})
