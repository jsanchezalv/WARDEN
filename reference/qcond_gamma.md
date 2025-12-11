# Conditional quantile function for gamma distribution

Conditional quantile function for gamma distribution

## Usage

``` r
qcond_gamma(rnd, shape, rate, lower_bound, s_obs)
```

## Arguments

- rnd:

  Vector of quantiles

- shape:

  The shape parameter

- rate:

  The rate parameter

- lower_bound:

  The lower bound to be used (current time)

- s_obs:

  is the survival observed up to lower_bound time, normally defined from
  time 0 as 1 - pgamma(q = lower_bound, rate, shape) but may be
  different if parametrization has changed previously

## Value

Estimate(s) from the conditional gamma distribution based on given
parameters

## Examples

``` r
qcond_gamma(rnd = 0.5, shape = 1.06178, rate = 0.01108,lower_bound = 1, s_obs=0.8)
#> [1] 87.94889
```
