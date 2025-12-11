# Draw from a gamma distribution based on mean and se

Draw from a gamma distribution based on mean and se

## Usage

``` r
rgamma_mse(n = 1, mean_v, se, seed = NULL)
```

## Arguments

- n:

  Number of draws (must be \>= 1)

- mean_v:

  A vector of the mean values

- se:

  A vector of the standard errors of the means

- seed:

  An integer which will be used to set the seed for this draw.

## Value

A single estimate from the gamma distribution based on given parameters

## Examples

``` r
rgamma_mse(n=1,mean_v=0.8,se=0.2)
#> [1] 1.050824
```
