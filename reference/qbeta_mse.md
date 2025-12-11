# Draw from a beta distribution based on mean and se (quantile)

Draw from a beta distribution based on mean and se (quantile)

## Usage

``` r
qbeta_mse(q, mean_v, se)
```

## Arguments

- q:

  Quantiles to be used

- mean_v:

  A vector of the mean values

- se:

  A vector of the standard errors of the means

## Value

A single estimate from the beta distribution based on given parameters

## Examples

``` r
qbeta_mse(q=0.5,mean_v=0.8,se=0.2)
#> [1] 0.8671142
```
