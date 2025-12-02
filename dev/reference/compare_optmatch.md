# Compares the equality of optmatch objects, ignoring attributes and group names.

This checks the equality of two optmatch objects. The only bits that
matter are unit names and the grouping. Other bits such as attributes,
group names, order, etc are ignored.

## Usage

``` r
compare_optmatch(o1, o2)
```

## Arguments

- o1:

  First optmatch object.

- o2:

  Second optmatch object.

## Value

TRUE if the two matches have the same memberships.

## Details

The names of the units can differ on any unmatched units, e.g., units
whose value in the optmatch object is `NA`. If matched objects have
differing names, this is automatically `FALSE`.

Note this ignores the names of the subgroups. So four members in
subgroups either `c("a", "a", "b", "b")` or `c("b", "b", "a", "a")`
would be identical to this call.
