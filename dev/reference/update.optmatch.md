# Performs an update on an `optmatch` object.

NB: THIS CODE IS CURRENTLY VERY MUCH ALPHA AND SOMEWHAT UNTESTED,
ESPECIALLY CALLING `update` ON AN OPTMATCH OBJECT CREATED WITHOUT AN
EXPLICIT DISTANCE ARGUMENT.

## Usage

``` r
# S3 method for class 'optmatch'
update(object, ...)
```

## Arguments

- object:

  `Optmatch` object to update.

- ...:

  Additional arguments to the call, or arguments with changed values.

## Value

An updated `optmatch` object.

## Details

Note that passing `data` again is strongly recommended. A warning will
be printed if the hash of the data used to generate the `optmatch`
object differs from the hash of the new `data`.

To obtain an updated call without performing the actual update, pass an
additional `evaluate = FALSE` argument.
