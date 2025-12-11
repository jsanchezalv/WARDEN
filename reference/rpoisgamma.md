# Draw time to event (tte) from a Poisson or Poisson-Gamma (PG) Mixture/Negative Binomial (NB) Process

Draw time to event (tte) from a Poisson or Poisson-Gamma (PG)
Mixture/Negative Binomial (NB) Process

## Usage

``` r
rpoisgamma(
  n,
  rate,
  theta = NULL,
  obs_time = 1,
  t_reps,
  seed = NULL,
  return_ind_rate = FALSE,
  return_df = FALSE
)
```

## Arguments

- n:

  The number of observations to be drawn

- rate:

  rate of the event (in terms of events per observation-time)

- theta:

  Optional. When omitted, the function simulates times for a Poisson
  process. Represents the shape of the gamma mixture distribution.
  Estimated and reported as theta in negative binomial regression
  analyses in r.

- obs_time:

  period over which events are observable

- t_reps:

  Optional. Number of TBEs to be generated to capture events within the
  observation window. When omitted, the function sets t_reps to the
  99.99th quantile of the Poisson (if no theta is provided) or negative
  binomial (if theta is provided). Thus, the risk of missing possible
  events in the observation window is 0.01%.

- seed:

  An integer which will be used to set the seed for this draw.

- return_ind_rate:

  A boolean that indicates whether an additional vector with the rate
  parameters used per observation is used. It will alter the structure
  of the results to two lists, one storing tte with name tte, and the
  other with name ind_rate

- return_df:

  A boolean that indicates whether a data.table object should be
  returned

## Value

Estimate(s) from the time to event based on poisson/Poisson-Gamma (PG)
Mixture/Negative Binomial (NB) distribution based on given parameters

## Details

Function to simulate event times from a Poisson or Poisson-Gamma (PG)
Mixture/Negative Binomial (NB) Process Event times are determined by
sampling times between events (TBEs) from an exponential distribution,
and cumulating these to derive the event times. Events occurring within
the set observation time window are retained and returned. For times for
a Poisson process, the provided rate is assumed constant. For a PG or
NB, the individual rates are sampled from a Gamma distribution with
shape = theta and scale = rate/theta.

## Examples

``` r
rpoisgamma(1,rate=1,obs_time=1,theta=1)
#> [[1]]
#> [1] 0.6935303
#> 
```
