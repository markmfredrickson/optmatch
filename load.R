## Use devtools to load most recent version
library(utils)
library("devtools")
devtools:::load_all()

## `utils` is loaded first to ensure proper ordering in the search
## stack, so that devtool's versions of `help` and `?` mask `utils`.
