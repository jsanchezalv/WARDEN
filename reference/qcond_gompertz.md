# Quantile function for conditional Gompertz distribution (lower bound only)

Quantile function for conditional Gompertz distribution (lower bound
only)

## Usage

``` r
qcond_gompertz(rnd, shape, rate, lower_bound = as.numeric(c(0)))
```

## Arguments

- rnd:

  Vector of quantiles

- shape:

  The shape parameter of the Gompertz distribution, defined as in the
  coef() output on a flexsurvreg object

- rate:

  The rate parameter of the Gompertz distribution, defined as in the
  coef() output on a flexsurvreg object

- lower_bound:

  The lower bound of the conditional distribution

## Value

Estimate(s) from the conditional Gompertz distribution based on given
parameters

## Examples

``` r
qcond_gompertz(rnd=0.5,shape=0.05,rate=0.01,lower_bound = 50)
#> [1] 5.007156
```
