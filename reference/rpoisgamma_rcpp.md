# Draw time to event (tte) from a Poisson or Poisson-Gamma (PG) Mixture/Negative Binomial (NB) Process using C++

Draw time to event (tte) from a Poisson or Poisson-Gamma (PG)
Mixture/Negative Binomial (NB) Process using C++

## Usage

``` r
rpoisgamma_rcpp(
  n,
  rate,
  theta = NULL,
  obs_time = 1,
  t_reps = NULL,
  seed = NULL,
  return_ind_rate = FALSE,
  return_df = FALSE
)
```

## Arguments

- n:

  The number of observations to be drawn

- rate:

  rate of the event (events per unit time)

- theta:

  Optional. If provided, Poisson-Gamma (NB). Represents gamma shape.

- obs_time:

  period over which events are observable

- t_reps:

  Optional. Number of TBEs to be generated to capture events within the
  observation window.

- seed:

  Optional integer seed for reproducibility.

- return_ind_rate:

  Logical: include individual rate vector in output when theta provided.

- return_df:

  Logical: return a data.frame with event-level rows (if TRUE).

## Value

If return_df=TRUE: a data.frame (or NULL if no events). Else: list with
tte and optionally ind_rate.

## Examples

``` r
rpoisgamma_rcpp(1, rate = 1, obs_time = 1, theta = 1)
#> $tte
#> $tte[[1]]
#> [1] 0.2308122 0.6734768 0.9150946
#> 
#> 


```
