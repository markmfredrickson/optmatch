require('optmatch')
data(nuclear, package="boot")
match_on(pr ~ cap, data = nuclear)
match_on(pr ~ date + cum.n, data = nuclear)
match_on(pr ~ date + cum.n, data = nuclear, within = exactMatch(pr ~ pt, data = nuclear))

if ( (require(splines)) )
  match_on(pr ~ ns(date,df=3) + cum.n, data = nuclear)

if (require(splines))
  match_on(pr ~ ns(date,df=3) + cum.n, data = nuclear, exactMatch(pr ~ pt, data = nuclear)
)
cum.n.q <- cut(nuclear$cum.n, quantile(nuclear$cum.n), include.lowest=TRUE)
match_on(pr ~ date + cum.n.q, data = nuclear)
match_on(pr ~ date + cum.n.q, data = nuclear, within = exactMatch(pr~pt, data = nuclear))

### should give error, incorrect mode
try(match_on(as.factor(pr) ~ cap, data = nuclear))
