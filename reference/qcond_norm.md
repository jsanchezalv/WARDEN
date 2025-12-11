# Conditional quantile function for normal distribution

Conditional quantile function for normal distribution

## Usage

``` r
qcond_norm(rnd, mean, sd, lower_bound, s_obs)
```

## Arguments

- rnd:

  Vector of quantiles

- mean:

  The mean parameter

- sd:

  The sd parameter

- lower_bound:

  The lower bound to be used (current time)

- s_obs:

  is the survival observed up to lower_bound time, normally defined from
  time 0 as 1 - pnorm(q = lower_bound, mean, sd) but may be different if
  parametrization has changed previously

## Value

Estimate(s) from the conditional normal distribution based on given
parameters

## Examples

``` r
qcond_norm(rnd = 0.5, mean = 1,sd = 1,lower_bound = 1, s_obs=0.8)
#> [1] 0.2533471
```
