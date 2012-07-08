require('optmatch')
data(nuclear, package="boot")
mahal.dist(pr~cap, nuclear)
mahal.dist(pr~cap, nuclear,
           inverse.cov=matrix(1,1,1,dimnames=list("cap", "cap")))
mahal.dist(pr~date+cum.n, nuclear)
mahal.dist(~date+cum.n, nuclear, pr~pt)
if ( (require(splines)) )
  mahal.dist(pr~ns(date,df=3)+cum.n, nuclear)
if (require(splines))
  mahal.dist(~ns(date,df=3)+cum.n, nuclear, pr ~ pt)
cum.n.q <- cut(nuclear$cum.n, quantile(nuclear$cum.n), include.lowest=TRUE)
mahal.dist(pr~date+cum.n.q, nuclear)
mahal.dist(~date+cum.n.q, nuclear, pr~pt)
### should give error, incorrect mode
try(mahal.dist(as.factor(pr)~cap, nuclear))

### should be OK:
mahal.dist(pr~date, nuclear, inverse.cov=diag(1))
### should complain about inverse.cov's lack of dimnames (but not give error)
mahal.dist(pr~date+cum.n, nuclear, inverse.cov=diag(2))
### same thing w/o complaint:
mahal.dist(pr~date+cum.n, nuclear, inverse.cov=structure(diag(2),
                                     dimnames=list(c("date","cum.n"),c("date","cum.n")))
           )
### also OK:
mahal.dist(pr~date+cum.n, nuclear, inverse.cov=structure(diag(2),
                                     dimnames=list(c("dateTRUE","cum.n"),c("dateTRUE","cum.n")))
           )

### but this should stop the show:
try(
    mahal.dist(pr~date+cum.n, nuclear, inverse.cov=structure(diag(2),
                                     dimnames=list(c("cum.n","date"),c("date","cum.n")))
               )
)
### and this should be another showstopper:
try(
    mahal.dist(pr~date+cum.n, nuclear, inverse.cov=structure(diag(2),
                                     dimnames=list(c("cum.n","date"),c("cum.n","date")))
               )
)
