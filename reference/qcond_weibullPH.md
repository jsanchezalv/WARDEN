# Conditional quantile function for WeibullPH (flexsurv)

Conditional quantile function for WeibullPH (flexsurv)

## Usage

``` r
qcond_weibullPH(rnd, shape, scale, lower_bound = as.numeric(c(0)))
```

## Arguments

- rnd:

  Vector of quantiles (between 0 and 1)

- shape:

  Shape parameter of WeibullPH

- scale:

  Scale (rate) parameter of WeibullPH (i.e., as in hazard = scale \*
  t^(shape - 1))

- lower_bound:

  Lower bound (current time)

## Value

Estimate(s) from the conditional weibullPH distribution based on given
parameters

## Examples

``` r
qcond_weibullPH(rnd = 0.5, shape = 2, scale = 0.01, lower_bound = 5)
#> [1] 4.711576
```
