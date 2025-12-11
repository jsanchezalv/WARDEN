# Survival Probaility function for conditional Gompertz distribution (lower bound only)

Survival Probaility function for conditional Gompertz distribution
(lower bound only)

## Usage

``` r
pcond_gompertz(time = 1, shape, rate, lower_bound = 0)
```

## Arguments

- time:

  Vector of times

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
pcond_gompertz(time=1,shape=0.05,rate=0.01,lower_bound = 50)
#> [1] 0.1174342
```
