context("MCFSolutions & co: S4 classes to encode min-cost-flow solutions")

test_that("Instantiation & validity", {
    expect_silent(spi1  <- new("SubProbInfo", data.frame(subproblem=c('a','b'), hashed_dist=c('a','b'),
                                           resolution=c(1,10), exceedance=c(.5, 2), CS_orig_dist=c(TRUE,FALSE),
                                           stringsAsFactors=F)
                               )
                  )
    expect_true(validObject(spi1))

    spi2  <- spi1
    colnames(spi2)[1]  <- "Subprob"
    expect_error(validObject(spi2), "Cols 1-5 should be")

    expect_silent(ni  <- new("NodeInfo",
                             data.frame(name='a', price=0.5, kind='bookkeeping',
                                        supply=1L, subproblem='b',
                                        stringsAsFactors=F)
                             )
                  )
    expect_error(new("NodeInfo",
                      data.frame(name='a', price=5L, kind='bookkeeping',
                                 supply=1L, subproblem='b',
                                 stringsAsFactors=F)
                     ),
                 "should have type double" # Not sure it's necessary, but insisting
                  )                        # that 'price' be double not integer

    expect_silent(ai  <- new("ArcInfo",
                             matches=data.frame(subproblem='a', treatment='b',
                                                control=c('c','d'),stringsAsFactors=F),
                             bookkeeping=data.frame(subproblem='a', startnode=c('c','d'),
                                                    endnode='(_Sink_)', flow=1L, stringsAsFactors=F)
                             )
                  )
    expect_error(new("ArcInfo",
                     matches=data.frame(subproblem='a', treatment='b',
                                        control=c('c','d'),stringsAsFactors=F),
                     bookkeeping=data.frame(subproblem='a', startnode=c('c','d'),
                                            endnode='(_Sink_)', flow=1.5, stringsAsFactors=F)
                     ), "should have type integer" # Not sure it's necessary, but insisting 
                  )                                # that 'flow' be integer not double
    
    expect_silent(mbls1  <- new("MatchablesInfo",
                                data.frame(name=c('a','b'), 'kind'=c('treatment', 'control'),
                                           subproblem=rep("c",2), stringsAsFactors=F)
                                )
                  )

    mbls2  <- mbls1
    colnames(mbls2)[3]  <- "Kind"
    expect_error(validObject(mbls2), "Columns should be", fixed=TRUE)

    mbls3  <- mbls1
    mbls3@.Data[[2]][1]  <- "Exposed"
    expect_error(validObject(mbls3), "values other than 'treatment' or 'control'", fixed=TRUE)

    expect_silent(new("MCFSolutions", subproblems=spi1,nodes=ni,arcs=ai,matchables=mbls1))
    expect_silent(mcf2  <- new("MCFSolutions", subproblems=spi1,nodes=ni,arcs=ai,matchables=mbls2))
    expect_error(validObject(mcf2, complete=TRUE), "Columns should be", fixed=TRUE)
})

test_that("c() methods", {

    spi1  <- new("SubProbInfo",
                 data.frame(subproblem=c('a','b'), hashed_dist=c('a','b'),
                            resolution=c(1,10), exceedance=c(.5, 2), CS_orig_dist=c(TRUE,FALSE),
                            stringsAsFactors=F)
                 )
    expect_silent(c(spi1, spi1))
    expect_silent(c(spi1, spi1, spi1))
    expect_silent(c(a=spi1, b=spi1)) # no confusion just b/c no `x=` arg!
    
    ni1  <- new("NodeInfo",
               data.frame(name='a', price=0.5, kind='bookkeeping',
                          supply=1L, subproblem='b',
                          stringsAsFactors=F)
               )
    expect_silent(c(ni1, ni1))
    expect_silent(c(ni1, ni1, ni1))
    
    mbls1  <- new("MatchablesInfo",
                  data.frame(name=c('a','b'), 'kind'=c('treatment', 'control'),
                             subproblem=rep("c",2), stringsAsFactors=F)
                  )
    expect_silent(c(mbls1, mbls1))

    ai  <- new("ArcInfo",
               matches=data.frame(subproblem='a', treatment='b',
                                  control=c('c','d'),stringsAsFactors=F),
               bookkeeping=data.frame(subproblem='a', startnode=c('c','d'),
                                      endnode='(_Sink_)', flow=1L, stringsAsFactors=F)
               )
    expect_silent(c(ai, ai))
    expect_silent(c(x=ai, y=ai, z=ai))

    mcf1  <- new("MCFSolutions", subproblems=spi1,nodes=ni1,arcs=ai,matchables=mbls1)
    expect_silent(c(mcf1, mcf1))
    expect_silent(c(y=mcf1, z=mcf1))
    expect_silent(c(mcf1, mcf1, mcf1))
})
