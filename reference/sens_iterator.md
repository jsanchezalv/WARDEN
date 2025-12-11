# Create an iterator based on sens of the current iteration within a scenario (DSA)

Create an iterator based on sens of the current iteration within a
scenario (DSA)

## Usage

``` r
sens_iterator(sens, n_sensitivity)
```

## Arguments

- sens:

  current analysis iterator

- n_sensitivity:

  total number of analyses to be run

## Value

Integer iterator based on the number of sensitivity analyses being run
and the total iterator

## Details

In a situation like a DSA, where two (low and high) scenarios are run,
sens will go from 1 to n_sensitivity\*2. However, this is not ideal as
the parameter selector may depend on knowing the parameter order (i.e.,
1, 2, 3...), which means resetting the counter back to 1 once sens
reaches n_sensitivity (or any multiple of n_sensitivity) is needed.

## Examples

``` r
sens_iterator(5,20)
#> [1] 5
sens_iterator(25,20)
#> [1] 5
```
