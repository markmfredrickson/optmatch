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

    propensity.model <- glm(z ~ w1 + w2, data = W, family =
      binomial())
    distances$propensity <- mdist(propensity.model)

The `glm` method is a wrapper around the `numeric` method for `mdist`. The
`numeric` method takes a vector of scores (for example, the linear prediction
for each unit from the model) and a vector indicating treatment status (`z`)
for each unit. This method returns the absolute difference between each
treated and control pair on their scores (additionally, 
the `glm` method rescales the data before invoking the `numeric` method). If
you wish to fit a "caliper" to your distance matrix, a hard limit on allowed
distances between treated and control units, you can pass a `caliper`
argument, a scalar numeric value. Any treated and control pair that is larger
than the caliper value will be replaced by `Inf`, an unmatchable value. The
`caliper` argument also applies to `glm` method. Calipers are covered in more
detail in the next sction.

The final convenience method of `mdist` is using an arbitrary function.  When
called, this function will receive  two `data.frame`s each of length `nt *
nc`, where `nt` is the number of treated units and `nc` is the number of
control units. The first argument represents the treated units, while the
second represents the control units. Combined the two variables represent every
possible treated-control pair. While this might sound complicated at first, it
is quite easy to use. For example, here is a function that sets distance to
infinity if the units differ on `w2` and the difference of `w1` otherwise:

    example.function <- function(t, c) {
      ifelse(t$w2 == c$w2, abs(t$w1 - c$w1), Inf)  
    }
    distances$fn <- mdist(example.function, data = W, z = W$z)

As in the previous examples, treated
and control pairs that have infinite distance will never be matched. 

### Combining and editing distances

We have created several representations of the matching problem, using
Euclidean distance, Mahalanobis distance, the estimated propensity score, and
an arbitrary function. We can combine these distances into single metric using
standard arithmetic functions:

    distances$all <- with(distances, euclid2 + mahal + propensity + fn)

You may find it convenient to work in smaller pieces at first and then stitch
the results together into a bigger distance. The `rbind` and `cbind` functions let us
add additional treated and control entries to a distance specification for
each of the existing control and treated units, respectively. For example, we
might want to combine a Mahalanobis score for units `n` through `s` with a
propensity score for units `t` through `z`:

    W.n.to.s <- W[c(LETTERS[1:13], letters[14:19]),]
    W.t.to.z <- W[c(LETTERS[1:13], letters[20:26]),]
    mahal.n.to.s <- mdist(z ~ w1 + w2, data = W)
    ps.t.to.z <- mdist(glm(z ~ w1 + w2, data = W.t.to.z, family = binomial()))
    distances$combined <- cbind(mahal.n.to.s, ps.t.to.z)

In a previous example, we used a function to compute conditional values:
infinity, if two units had different values of `w2`; the difference of `w1`
otherwise. `optmatch` has several functions that follow this basic pattern:
infinity, under some condition; zero otherwise. Since distances can be added,
we can recreate the results of the function based `mdist` call:

    tmp <- sqrt(mdist(z ~ w1, data = W, inv.scale.matrix = 1))
    distances$exact <- tmp + exactMatch(z ~ w2, data = W)
    all(distances$fn == as.matrix(distances$exact))

The `exactMatch` function creates "stratified" matching problems, in which
there are subgroups that are completely separate. Such matching problems are
often much easier to solve than problems where a treated unit could be
connected to any control unit.

The use of `as.matrix` in the comparison is because the results of
`exactMatch` are a special kind of matrix that only records the finite entries
and ignores any entries that would otherwise be infinity. In our small
example, this does not waste much time or space, but in bigger problems, where many
treatment-control pairs are disallowed, this data type can be very beneficial.

There is another method for creating reduced matching problems. The `caliper`
function compares each entry in an existing distance specification and
disallows any that are larger than a specified value. For example, we can trim
our previous combined distance to anything smaller than the median value:

    distances$median.caliper <- caliper(distances$all, median(distances$all))
    distances$all.trimmed <- with(distances, all + median.caliper)

Like the `exactMatch` function, the results of `caliper` used the sparse
matrix representation mentioned above, so can be very efficient for large,
sparse problems. As noted previously, if using the `glm` or `numeric` methods
of `mdist`, passing the caliper's width in the `caliper` argument can be more
efficient.

### Speeding up computation

In addition to the space advantages of only storing the finite entries in a
sparse matrix, the results of `exactMatch` and `caliper` can be used to speed
up computation of *new* distances. The `mdist` function that we saw earlier
has an argument called `within` that helps filter the resulting
computation to only the finite entries in the `within` matrix. Since `exactMatch` and `caliper`
use finite entries denote valid pairs, they make excellent sources of
the `within` argument.

Instead of creating the entire Euclidean distance matrix and *then* filtering
out cross-strata matches, we use the results of `exactMatch` to compute only
the interesting cases:

    tmp <- exactMatch(z ~ w2, data = W) # same as before
    distances$exact2 <- sqrt(mdist(z ~ w1, data = W, inv.scale.matrix = 1,
    within = tmp))
    identical(distances$exact2, distances$exact)

Users of previous versions of `optmatch` may notice that the `within`
argument is similar to the old `structure.formula` argument. Like
`within`, `structure.formula` focused distance on within strata pairs.
Unlike `structure.formula`, the `within` argument allows using *any*
distance specification as an argument, including those created with `caliper`. For
example, here is the Mahalanobis distance computed only for units that differ
by less than one on the propensity score.

    distances$mahal.trimmed <- mdist(z ~ w1 + w2, data = W,
      exclusions = mdist(propensity.model, caliper = 1))

### Generating the match

Now that we have generated several distances specifications, let's put them to
use. Here is the simplest way to evaluate all distances specifications:

    matches <- lapply(distances, function(x) { fullmatch(x, data = W) })

The result of the matching process is a named factor, where the names
correspond to the units (both treated and control) and the levels of the
factors are the matched groups. Including the `data` argument is highly
recommended. This argument will make sure that the result of `fullmatch` will
be in the same order as the original `data.frame` that was used to build the
distance specification. This will make appending the results of `fullmatch`
on to the original `data.frame` much more convenient.

The `fullmatch` function as several arguments for fine tuning the allowed
ratio of treatment to control units in a match, and how much of the pool to
throw away as unmatchable. One common pattern for these arguments are pairs:
one treated to one control unit. Not every distance specification is amendable
to this pattern (e.g. when there are more treated units than control units in
`exactMatch` created stratum). However, it can be done with the Mahalanobis
distance matrix we created earlier:

    mahal.match <- pairmatch(distances$mahal, data = W)

Like `fullmatch`, `pairmatch` also allows fine tuning the ratio of matches to
allow larger groupings. It is can be helpful as it computes what percentage of
the group to throw away, giving better odds of successfully finding a matching
solution.

Once one has generated a match, you may wish to view the results. The results
of calls to `fullmatch` or `pairmatch` produce `optmatch` objects (specialized
factors). This object has a special option to the `print` method which groups
the units by factor level:

    print(mahal.match, grouped = T)

If you wish to join the match factor back to the original `data.frame`: 

    W.matched <- cbind(W, matches = mahal.match)

Make sure to include the `data` argument to `fullmatch` or `pairmatch`,
otherwise results are not guaranteed to be in the same order as your original
`data.frame` or `matrix`. 


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
[Windows](http://cran.r-project.org/bin/windows/Rtools/). You will also need a
copy of GNU `make` to create the package from source (standard on Linux,
included with [Apple's developer tools](http://developer.apple.com), included
with the [Cygwin](http://www.cygwin.com/) UNIX tools for Windows).

In `R`, you will need to install and load the `roxygen2` package:

    > install.packages("roxygen2")

You will need to download the branch you wish to install. You can navigate to
the [Optmatch project page](http://github.com/markmfredrickson/optmatch) to
find the branch you wish to use. Switch to the proper branch and use the zip
download button. Unzip the package. In a terminal window:

    $ cd /path/to/package
    $ make package

This should build a `optmatch_VERSION.tar.gz` file. You can install it in a
local directory (for example `~/R/optmatch.demo`) using:

    $ mkdir -p ~/R/optmatch.demo
    $ R CMD Install --no-multiarch --library=~/R/optmatch.demo ./optmatch_VERSION.tar.gz

You can then load the library in `R` using:

    > library("optmatch", lib.loc = "~/R/optmatch.demo") 

### Developing for Optmatch

We welcome patches to add features, fix bugs, or otherwise improve the package.
To develop `optmatch`, you will need to have a working installation of `git`
and all the software mentioned in the previous section. Instead of downloading
the source directly, fork the project and github and clone a working copy from
your forked project:

    $ git clone git@github.com:YOURUSERNAME/optmatch.git

We prefer changes that include unit tests demonstrating the problem or showing
how the new feature should be added. The test suite uses the
[testthat](http://github.com/hadley/test_that) package to write and run tests.
See the `inst/tests` directory for examples. To run the test suite, use:

    $ make test

New features should include documentation. We prefer inline
[Roxygen](http://roxygen.org/) style documentation to raw `.Rd` files. Any
`Makefile` task that builds the package will create the documentation. These
tasks are also useful for development:

- `R`: starts up an interactive session with `optmatch` loaded.
- `spell`: checks the spelling of all Rd files in the package using the
  []`aspell_*`
  functions](http://stat.ethz.ch/R-manual/R-patched/library/utils/html/aspell-utils.html)
  in the 'utils' package. - `check`: runs `R CMD Check` on a built package, includes the `test` and `spell` tasks
- `package`: creates a `tar.gz` of the package
- `check`: runs `R CMD check` on the package

If you need to edit the `DESCRIPTION` file (e.g. to add a suggested package),
edit `DESCRIPTION.template` instead. This file is used during the build process
to insert the version and date information.

You will need to have the `aspell` program installed and available on your
`$PATH` in order to use the `spell` task. If you introduce new words that are
not coveraged in the standard English dictionary, you can add them to the
`lexicon.txt` file. There is also a task to create the lexicon from all found
misspelled words (`make lexicon.txt`).  For more information on the use of
spell checkers in R, see [Watch your
spelling!](http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Hornik+Murdoch.pdf)
in the R journal.

When your change is ready, make a pull request on github.

