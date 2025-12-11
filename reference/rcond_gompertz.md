# Draw from a conditional Gompertz distribution (lower bound only)

Draw from a conditional Gompertz distribution (lower bound only)

## Usage

``` r
rcond_gompertz(n = 1, shape, rate, lower_bound = 0, seed = NULL)
```

## Arguments

- n:

  The number of observations to be drawn

- shape:

  The shape parameter of the Gompertz distribution, defined as in the
  coef() output on a flexsurvreg object

- rate:

  The rate parameter of the Gompertz distribution, defined as in the
  coef() output on a flexsurvreg object

- lower_bound:

  The lower bound of the conditional distribution

- seed:

  An integer which will be used to set the seed for this draw.

## Value

Estimate(s) from the conditional Gompertz distribution based on given
parameters

## Examples

``` r
rcond_gompertz(1,shape=0.05,rate=0.01,lower_bound = 50)
#> [1] 1.42697
```
