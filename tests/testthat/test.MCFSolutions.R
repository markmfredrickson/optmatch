context("MCFSolutions & co: S4 classes to encode min-cost-flow solutions")

test_that("Instantiation & validity", {
    expect_silent(spi1  <- new("SubProbInfo", data.frame(groups=c('a','b'), flipped=logical(2), hashed_dist=c('a','b'),
                                           resolution=c(1,10), exceedance=c(.5, 2), CS_orig_dist=c(TRUE,FALSE),
                                           stringsAsFactors=F)
                               )
                  )
    expect_true(validObject(spi1))

    spi2  <- spi1
    colnames(spi2)[1]  <- "Subprob"
    expect_error(validObject(spi2), "Cols 1-6 should be")

    expect_silent(ni1  <- new("NodeInfo",
                             data.frame(name='a', price=0.5, upstream_not_down=TRUE,
                                        supply=1L, groups = as.factor('b'),
                                        stringsAsFactors=F)
                             )
                  )
    expect_silent(ni1f  <- new("NodeInfo",
                              data.frame(name=c('a', '(_Sink_)', '(_End_)'), price=0.5,
                                         upstream_not_down=c(FALSE, NA, NA),
                                         supply=c(0L,-1L,-2L), groups = as.factor('b'),
                                         stringsAsFactors=F)
                              )
                  )
    expect_error(new("NodeInfo",
                      data.frame(name='a', price='foo', upstream_not_down=FALSE,
                                 supply=1L, groups = as.factor('b'),
                                 stringsAsFactors=F)
                     ),
                 "should be a numeric" # Not sure it's necessary, but insisting
                  )                        # that 'price' be double not integer
    expect_silent(ai  <- new("ArcInfo",
                             matches=data.frame(groups = as.factor('a'), upstream = as.factor('b'),
                                                downstream = as.factor(c('c','d')), stringsAsFactors=F),
                             bookkeeping=data.frame(groups = as.factor('a'), start = as.factor(c('c','d')),
                                                    end = as.factor('(_Sink_)'), flow=1L,
                                                    capacity=1L, stringsAsFactors=F)
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
                     matches=data.frame(groups = as.factor('a'), upstream = as.factor('b'),
                                        downstream=c('c','d'),stringsAsFactors=F),
                     bookkeeping=data.frame(groups='a', start=c('c','d'),
                                            end='(_Sink_)', flow=2L,
                                            capacity=1L, stringsAsFactors=F)
                     ), "flow can be now greater than capacity" 
                  )                                
    
    expect_silent(mbls1  <- new("MatchablesInfo",
                                data.frame(name=c('a','b'), 'row_unit'=c(TRUE, FALSE),
                                           groups = as.factor(rep("a",2)), stringsAsFactors=F)
                                )
                  )

    mbls2  <- mbls1
    colnames(mbls2)[3]  <- "Kind"
    expect_error(validObject(mbls2), "Columns should be", fixed=TRUE)    
    
    expect_silent(mcf1  <- new("MCFSolutions", subproblems=spi1,nodes=ni1,arcs=ai,matchables=mbls1))
    expect_silent(mcf2  <- new("MCFSolutions", subproblems=spi1,nodes=ni1,arcs=ai,matchables=mbls2))
    expect_error(validObject(mcf2, complete=TRUE), "Columns should be", fixed=TRUE)
    mbls3  <- new("MatchablesInfo",
                  data.frame(name=c('a','b'), 'row_unit'=c(TRUE, FALSE),
                             groups = as.factor(rep("d",2)), stringsAsFactors=F)
                  )
    mcf3  <- mcf1
    mcf3@matchables  <- mbls3 #mismatch between subproblems here vs elswhere in object
    expect_error(validObject(mcf3), "etected subproblems")

    expect_silent(mcf1f  <- new("MCFSolutions", subproblems=spi1,nodes=ni1f,arcs=ai,matchables=mbls1))
    expect_silent(as(mcf1f, "FullmatchMCFSolutions"))
})

test_that("c() methods", {
    spi1  <- new("SubProbInfo",
                 data.frame(groups=c('a','b'), flipped=logical(2), hashed_dist=c('a','b'),
                            resolution=c(1,10), exceedance=c(.5, 2), CS_orig_dist=c(TRUE,FALSE),
                            stringsAsFactors=F)
                 )
    spi2  <- new("SubProbInfo",
                 data.frame(groups=c('c'), flipped=logical(1), hashed_dist=c('a'),
                            resolution=c(1), exceedance=c(.5), CS_orig_dist=c(TRUE),
                            stringsAsFactors=F)
                 )
    spi3  <- new("SubProbInfo",
                 data.frame(groups=c('d'), flipped=logical(1), hashed_dist=c('a'),
                            resolution=c(1), exceedance=c(.5), CS_orig_dist=c(TRUE),
                            stringsAsFactors=F)
                 )
    
    expect_silent(c(spi1, spi2))
    expect_silent(c(spi1, spi2, spi3))
    expect_silent(c(a=spi1, b=spi2)) # no confusion just b/c no `x=` arg!
    
    ni1  <- new("NodeInfo",
               data.frame(name='a', price=0.5, upstream_not_down=TRUE,
                          supply=1L, groups= as.factor('b'),
                          stringsAsFactors=F)
               )
    expect_silent(c(ni1, ni1))
    expect_silent(c(ni1, ni1, ni1))
    ni2  <- new("NodeInfo",
                data.frame(name=c('a', '(_Sink_)', '(_End_)'), price=0.5,
                           upstream_not_down=c(FALSE, NA, NA),
                          supply=c(0L,-1L,-2L), groups=as.factor('c'),
                          stringsAsFactors=F)
               )
    ni1ni2 <- c(ni1, ni2)
    expect_equal(ni1ni2$name, c("a", "a", "(_Sink_)", "(_End_)"))
    expect_equal(levels(ni1ni2$groups), c("b", "c"))

    mbls1  <- new("MatchablesInfo",
                  data.frame(name=c('a','b'), 'row_unit'=c(TRUE, FALSE),
                             groups = as.factor(rep("a",2)), stringsAsFactors=F)
                  )
    expect_silent(c(mbls1, mbls1))
    mbls2  <- new("MatchablesInfo",
                  data.frame(name=c('a','b'), 'row_unit'=c(TRUE, FALSE),
                             groups = as.factor(rep("c",2)), stringsAsFactors=F)
                  )

    mb1mb2 <- c(mbls1, mbls2)
    expect_equal(mb1mb2$name,  c('a', 'b', 'a', 'b'))
    expect_equal(levels(mb1mb2$groups), c('a', 'c'))


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


    mcf1  <- new("MCFSolutions", subproblems=spi1, nodes=ni1, arcs=ai1, matchables=mbls1)
    expect_error(c(mcf1, mcf1), "uplicates")


    mcf2 <- new("MCFSolutions", subproblems=spi2,nodes=ni2,arcs=ai2,matchables=mbls2)
    
    expect_silent(c(mcf1, mcf2))
    expect_silent(c(y=mcf1, z=mcf2))

    mcf3  <- mcf2
    mcf3@matchables  <-
        new("MatchablesInfo",
                  data.frame(name=c('a','b'), 'row_unit'=c(TRUE, FALSE),
                             groups = as.factor(rep("d",2)), stringsAsFactors=F)
            )
    ## this `mcf3` is now corrupted:
    expect_error(validObject(mcf3))
    ## however, to save time c() just assumes each 
    ## of the MCFSolutions objects it's combining is individually valid. 
    expect_silent(c(mcf1, mcf3))

    mcf2f  <- as(mcf2, "FullmatchMCFSolutions")
    expect_is(c(mcf2f, mcf1), "FullmatchMCFSolutions")
    expect_is(c(mcf1, mcf2f), "MCFSolutions")
})
