# Conditional quantile function for weibull distribution

Conditional quantile function for weibull distribution

## Usage

``` r
qcond_weibull(rnd, shape, scale, lower_bound = as.numeric(c(0)))
```

## Arguments

- rnd:

  Vector of quantiles

- shape:

  The shape parameter as in R stats package weibull

- scale:

  The scale parameter as in R stats package weibull

- lower_bound:

  The lower bound to be used (current time)

## Value

Estimate(s) from the conditional weibull distribution based on given
parameters

## Examples

``` r
qcond_weibull(rnd = 0.5,shape = 3,scale = 66.66,lower_bound = 50)
#> [1] 19.12624
```
