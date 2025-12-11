# Conditional quantile function for lognormal distribution

Conditional quantile function for lognormal distribution

## Usage

``` r
qcond_lnorm(rnd, meanlog, sdlog, lower_bound, s_obs)
```

## Arguments

- rnd:

  Vector of quantiles

- meanlog:

  The meanlog parameter

- sdlog:

  The sdlog parameter

- lower_bound:

  The lower bound to be used (current time)

- s_obs:

  is the survival observed up to lower_bound time, normally defined from
  time 0 as 1 - plnorm(q = lower_bound, meanlog, sdlog) but may be
  different if parametrization has changed previously

## Value

Estimate(s) from the conditional lognormal distribution based on given
parameters

## Examples

``` r
qcond_lnorm(rnd = 0.5, meanlog = 1,sdlog = 1,lower_bound = 1, s_obs=0.8)
#> [1] 2.502045
```
