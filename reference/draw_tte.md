# Draw a time to event from a list of parametric survival functions

Draw a time to event from a list of parametric survival functions

## Usage

``` r
draw_tte(
  n_chosen,
  dist,
  coef1 = NULL,
  coef2 = NULL,
  coef3 = NULL,
  ...,
  beta_tx = 1,
  seed = NULL
)
```

## Arguments

- n_chosen:

  The number of observations to be drawn

- dist:

  The distribution; takes values
  'lnorm','norm','mvnorm','weibullPH','weibull','llogis','gompertz','gengamma','gamma','exp','beta','poisgamma'

- coef1:

  First coefficient of the distribution, defined as in the coef() output
  on a flexsurvreg object (rate in "rpoisgamma")

- coef2:

  Second coefficient of the distribution, defined as in the coef()
  output on a flexsurvreg object (theta in "rpoisgamma")

- coef3:

  Third coefficient of the distribution, defined as in the coef() output
  on a flexsurvreg object (not used in "rpoisgamma")

- ...:

  Additional arguments to be used by the specific distribution (e.g.,
  return_ind_rate if dist = "poisgamma")

- beta_tx:

  Parameter in natural scale applied in addition to the scale/rate
  coefficient -e.g., a HR if used in an exponential- (not used in
  "rpoisgamma" nor "beta")

- seed:

  An integer which will be used to set the seed for this draw.

## Value

A vector of time to event estimates from the given parameters

## Details

Other arguments relevant to each function can be called directly

## Examples

``` r
draw_tte(n_chosen=1,dist='exp',coef1=1,beta_tx=1)
#> [1] 0.329793
draw_tte(n_chosen=10,"poisgamma",coef1=1,coef2=1,obs_time=1,return_ind_rate=FALSE)
#> [[1]]
#> numeric(0)
#> 
#> [[2]]
#> [1] 0.6901357
#> 
#> [[3]]
#> numeric(0)
#> 
#> [[4]]
#> [1] 0.7689012
#> 
#> [[5]]
#> numeric(0)
#> 
#> [[6]]
#> [1] 0.4775767 0.5404240 0.5600788 0.6349382 0.6412048 0.6822714
#> 
#> [[7]]
#> numeric(0)
#> 
#> [[8]]
#> numeric(0)
#> 
#> [[9]]
#> [1] 0.1277729
#> 
#> [[10]]
#> [1] 0.3075535 0.9511123
#> 
```
