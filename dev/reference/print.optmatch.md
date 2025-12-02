# Printing `optmatch` objects.

Printing `optmatch` objects.

## Usage

``` r
# S3 method for class 'optmatch'
print(x, quote = FALSE, grouped = FALSE, ...)
```

## Arguments

- x:

  The `optmatch` object, as returned by
  [`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md)
  or
  [`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md).

- quote:

  A boolean indicating if the matched group names should be quoted or
  not (default is not to quote).

- grouped:

  A logical indicating if the object should printed in the style of a
  named `factor` object (`grouped = TRUE`) or as a table of group names
  and members.

- ...:

  Arguments passed to
  [`print.default`](https://rdrr.io/r/base/print.default.html).

## See also

[`fullmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/fullmatch.md),
[`pairmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/pairmatch.md),
[`print`](https://rdrr.io/pkg/RItools/man/print.xbal.html),
[`summary.optmatch`](https://markmfredrickson.github.io/optmatch/dev/reference/optmatch.md)

## Examples

``` r
data(nuclearplants)
fm <- fullmatch(pr ~ t1 + t2, data = nuclearplants)

print(fm)
#>    H    I    A    J    B    K    L    M    C    N    O    P    Q    R    S    T 
#>  1.3  1.1  1.1  1.8  1.2  1.3  1.3  1.3  1.3  1.8  1.8  1.3  1.5  1.3  1.9 1.10 
#>    U    D    V    E    W    F    X    G    Y    Z    d    e    f    a    b    c 
#>  1.7  1.4  1.4  1.5  1.2  1.6  1.5  1.7  1.3  1.6  1.3  1.3  1.3  1.8  1.9 1.10 
print(fm, grouped = TRUE)
#>  Group                         Members
#>    1.1                            I, A
#>   1.10                            T, c
#>    1.2                            B, W
#>    1.3 H, K, L, M, C, P, R, Y, d, e, f
#>    1.4                            D, V
#>    1.5                         Q, E, X
#>    1.6                            F, Z
#>    1.7                            U, G
#>    1.8                      J, N, O, a
#>    1.9                            S, b
```
