# Optmatch: Optimal Fullmatching for R

The `optmatch` package implements the optimal fullmatching algorithm for
bipartite matching problems. Given a matrix describing the distances between
two groups (where one group is represented by row entries, and the other by
column entries), the algorithm finds a matching between units that
minimizinges the average within groupd distances. This algorithm is a popular
choice for covariate balancing applications (e.g. propensity score matching),
but it also can be useful for design stage applications such as blocking. For
more on the application and its implementation, see:

    Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and
     related designs via network flows, JCGS 15 609-627.

`optmatch` is available on [CRAN](http://cran.r-project.org):

    > install.packages("optmatch")
    > library("optmatch")

##  Using a development version of Optmatch

This section will help you get the latest development version of `optmatch` and
start using the latest features. Before starting, you should know which branch
you wish to install. Currently, the "master" branch is the main code base. The 
"s4" branch contains a significant update to the class hierachy to use `R`'s 
so-called "S4" classes.

### Fetching and installing in a local directory

Chances are you already have an installation of `optmatch` that you use. These
directions will install the development version in a way that will not
overwrite your existing installation.

You must have the Fortran extensions for package building included. These can
be had from CRAN: [OS X](http://cran.r-project.org/bin/macosx/tools/),
[Windows](http://cran.r-project.org/bin/windows/Rtools/).

In `R`, you will need to install and load the `devtools` package:

    > install.packages("devtools")
    > library("devtools")

Next, pick a location to install the package. For example, I created a
directory called `~/R/optmatch.demo/` (`~` is short for my home directory on a
UNIX system). For this session, we will set the library path to look in this
location first and install the package there:

    > .libPaths("~/R/optmatch.demo/") # <- your path here
    > install_github("optmatch", user = "markmfredrickson", branch = "s4")

The function `install_github` will load the package automatically. In the
future, if you wish load the downloaded version of `optmatch` in a new `R`
session you can use this one-liner:

    > library("optmatch", lib.loc = "~/R/optmatch.demo") # <- your path here



