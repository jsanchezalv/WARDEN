# Conditional quantile function for exponential distribution

Conditional quantile function for exponential distribution

## Usage

``` r
qcond_exp(rnd, rate)
```

## Arguments

- rnd:

  Vector of quantiles

- rate:

  The rate parameter

  Note taht the conditional quantile for an exponential is independent
  of time due to constant hazard

## Value

Estimate(s) from the conditional exponential distribution based on given
parameters

## Examples

``` r
qcond_exp(rnd = 0.5,rate = 3)
#> [1] 0.2310491
```
