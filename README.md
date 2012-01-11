# Optmatch: Optimal Fullmatching for R

The `optmatch` package implements the optimal full matching algorithm for
bipartite matching problems. Given a matrix describing the distances between
two groups (where one group is represented by row entries, and the other by
column entries), the algorithm finds a matching between units that
minimizes the average within grouped distances. This algorithm is a popular
choice for covariate balancing applications (e.g. propensity score matching),
but it also can be useful for design stage applications such as blocking. For
more on the application and its implementation, see:

    Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and
     related designs via network flows, JCGS 15 609-627.

`optmatch` is available on [CRAN](http://cran.r-project.org):

    > install.packages("optmatch")
    > library("optmatch")

## Using Optmatch

In addition to the optimal full matching algorithm, the package contains
useful functions for generating distance specifications, combining and editing
distance specifications, and summarizing and displaying matches. This walk
through shows how to use these tools in your matching workflow.

### Simulated data

Before we start, let's generate some simulated data. We will have two groups,
the "treated" and "control" groups. Without our knowledge, nature assigned
units from a pool into one of these two groups. The probability of being a
treated unit depends on some covariates. In the vector `Z`, let a 1 denote
treated units and 0 denote control units

    set.seed(20120111) # set this to get the exact same answers as I do
    n <- 26 # chosen so we can divide the alphabet in half
    W <- data.frame(w1 = rbeta(n, 4, 2), w2 = rbinom(n, 1, p = .33))

    # nature assigns to treatment
    tmp <- numeric(n)
    tmp[sample(1:n, prob = W$w1^(1 + W$w2), size = n/2)] <- 1
    W$z <- tmp

    # for convenience, let's give the treated units capital letter names
    tmp <- character(n)
    tmp[W$z == 1] <- LETTERS[1:(n/2)]
    tmp[W$z == 0] <- letters[(26 - n/2 + 1):26]
    rownames(W) <- tmp

As we can see with a simple table and plot, these groups are not balanced on
the covariates, as they would be (in expectation) with a randomly assigned
treatment.

    table(W$w2, W$z)
    library(lattice) ; densityplot(W$w1, groups = W$z)

The next steps use the covariates to pair up similar treated and control
units. For more on assessing the amount and severity of imbalance between
groups on observed covariates, see the
[RItools](http://github.com/markmfredrickson/RItools) `R` package.

### Setting up distances

These two groups are different, but how different are individual treated units
from individual control units? In answering this question, we will produce
several distance specifications: matrices of treated units (rows) by control
units (columns) with entries denoting distances. `optmatch` provides several
ways of generating these matrices so that you don't have to do it by hand.

Let's begin with a simple (squared) Euclidean distance on the space defined by `W`:

    distances <- list()
    distances$euclid2 <- mdist(z ~ w1 + w2, data = W, inv.scale.matrix =
      diag(2))

The `inv.scale.matrix` argument tells the mdist function how to adjust the
scales of the different covariates. In the case of Euclidean distance, we
don't want any scaling, so we use a 2 by 2 identity matrix. The result of this
computation is a matrix, and we can get the non-squared Euclidean distance
using the usual function:

    class(distances$euclid2)
    distances$euclid <- sqrt(distances$euclid2)

While the scales of a uniform and a Bernoulli variable are similar, you may
find that you have different variables that are on wildly different scales or
that covary. A better choice might be to use the so called (squared)
[Mahalanobis distance](http://en.wikipedia.org/wiki/Mahalanobis_distance):

    distances$mahal <- mdist(z ~ w1 + w2, data = W)

The default argument for `inv.scale.matrix` is the inverse covariance matrix
of the data. The distance is calculated as:

  (T - C)^T S^(-1) (T - C)^T

So you can see how other scale matrices can be used.

To create distances, we could also try regressing the treatment indicator on
the covariates and computing the difference distance for each treated and
control pair. To make this process easier, `mdist` has methods for `glm`
objects (and for big data problems, `bigglm` objects):

    distances$propensity <- mdist(glm(z ~ w1 + w2, data = W, family =
      binomial()))

The final convenience method of `mdist` is using an arbitrary function.  When
called, this function will receive a two `data.frame`s each of length `nt *
nc`, where `nt` is the number of treated units and `nc` is the number of
control units. The first argument represents the treated units, while the
second represents the control units. Combined the two vectors represents every
possible treated-control pair. While this might sound complicated at first, it
is quite easy to use. For example, here is a function that sets distance to
infinity if the units differ on `w2` and the difference of `w1` otherwise:

    example.function <- function(t, c) {
      ifelse(t$w2 == c$w2, abs(t$w1 - c$w1), Inf)  
    }
    distances$w1diff <- mdist(example.function, data = W, z = W$z)

This final example shows useful property of distance specification: treated
and control pairs that have infinite distance will never be matched! We will
show how to use this property to fine tune matching problems in the next
section.

### Combining and editing distances

### Speeding up computation

### Generating the match

### Evaluating candidate matches


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



