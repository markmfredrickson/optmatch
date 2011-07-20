require('optmatch')
data(nuclear, package="boot")
mdist(pr ~ cap, data = nuclear)
mdist(pr ~ cap, data = nuclear,
           inv.scale.matrix = matrix(1,1,1, dimnames=list("cap", "cap")))
mdist(pr ~ date + cum.n, data = nuclear)
mdist(pr ~ date + cum.n, data = nuclear, exclusions = exactMatch(pr ~ pt, data = nuclear))

if ( (require(splines)) )
  mdist(pr ~ ns(date,df=3) + cum.n, data = nuclear)

if (require(splines))
  mdist(pr ~ ns(date,df=3) + cum.n, data = nuclear, exactMatch(pr ~ pt, data = nuclear)
)
cum.n.q <- cut(nuclear$cum.n, quantile(nuclear$cum.n), include.lowest=TRUE)
mdist(pr ~ date + cum.n.q, data = nuclear)
mdist(pr ~ date + cum.n.q, data = nuclear, exclusions = exactMatch(pr~pt, data = nuclear))

### should give error, incorrect mode
try(mdist(as.factor(pr) ~ cap, data = nuclear))
