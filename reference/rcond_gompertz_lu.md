# Draw from a Conditional Gompertz distribution (lower and upper bound)

Draw from a Conditional Gompertz distribution (lower and upper bound)

## Usage

``` r
rcond_gompertz_lu(
  n,
  shape,
  rate,
  lower_bound = 0,
  upper_bound = Inf,
  seed = NULL
)
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

- upper_bound:

  The upper bound of the conditional distribution

- seed:

  An integer which will be used to set the seed for this draw.

## Value

Estimate(s) from the Conditional Gompertz distribution based on given
parameters

## Examples

``` r
rcond_gompertz_lu(1,shape=0.05,rate=0.01,lower_bound = 50)
#> [1] 10.10558
```
