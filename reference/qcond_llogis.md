# Conditional quantile function for loglogistic distribution

Conditional quantile function for loglogistic distribution

## Usage

``` r
qcond_llogis(rnd, shape, scale, lower_bound = as.numeric(c(0)))
```

## Arguments

- rnd:

  Vector of quantiles

- shape:

  The shape parameter

- scale:

  The scale parameter

- lower_bound:

  The lower bound to be used (current time)

## Value

Estimate(s) from the conditional loglogistic distribution based on given
parameters

## Examples

``` r
qcond_llogis(rnd = 0.5,shape = 1,scale = 1,lower_bound = 1)
#> [1] 2
```
