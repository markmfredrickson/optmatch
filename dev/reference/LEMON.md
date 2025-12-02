# (Internal) Helper function for accessing algorithms in LEMON solver

(Internal) Helper function for accessing algorithms in LEMON solver

## Usage

``` r
LEMON(algorithm = "CycleCancelling")
```

## Arguments

- algorithm:

  LEMON algorithm to use. Choices are "CycleCancelling",
  "CapacityScaling", "CostScaling", "NetworkSimplex". Default is
  "CycleCancelling".

## Value

String of the form "LEMON.\<algorithm\>"
