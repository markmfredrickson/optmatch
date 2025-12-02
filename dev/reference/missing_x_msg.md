# (Internal) If the x argument does not exist for match_on, fullmatch, or pairmatch, use this function to print a helpful message.

(Internal) If the x argument does not exist for match_on, fullmatch, or
pairmatch, use this function to print a helpful message.

## Usage

``` r
missing_x_msg(x_str, data_str, ...)
```

## Arguments

- x_str:

  x as a string of code, usually deparse(substitute(x))

- data_str:

  data arg as string of code

- ...:

  will look for 'z = \<stuff\>' in the extra args of caller

## Value

string a helpful error message

## Author

Josh Buckner
