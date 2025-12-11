# Creates a vector of indicators (0 and 1) for sensitivity/DSA analysis

Creates a vector of indicators (0 and 1) for sensitivity/DSA analysis

## Usage

``` r
create_indicators(sens, n_sensitivity, elem, n_elem_before = 0)
```

## Arguments

- sens:

  current analysis iterator

- n_sensitivity:

  total number of analyses to be run

- elem:

  vector of 0s and 1s of elements to iterate through (1 = parameter is
  to be included in scenario/DSA)

- n_elem_before:

  Sum of 1s (# of parameters to be included in scenario/DSA) that go
  before elem

## Value

Numeric vector composed of 0 and 1, where value 1 will be used by
`pick_val_v` to pick the corresponding index in its `sens` argument

## Details

n_elem_before is to be used when several indicators want to be used
(e.g., for patient level and common level inputs) while facilitating
readibility of the code

## Examples

``` r
create_indicators(10,20,c(1,1,1,1))
#> [1] 0 0 0 0
create_indicators(7,20,c(1,0,0,1,1,1,0,0,1,1),2)
#>  [1] 0 0 0 0 0 0 0 0 1 0
```
