# Draw Time-to-Event with Time-Dependent Covariates and Luck Adjustment

Simulate a time-to-event (TTE) from a parametric distribution with
parameters varying over time. User provides parameter functions and
distribution name. The function uses internal survival and conditional
quantile functions, plus luck adjustment to simulate the event time. See
the vignette on avoiding cycles for an example in a model.

## Usage

``` r
qtimecov(
  luck,
  a_fun,
  b_fun = NULL,
  dist = "exp",
  dt = 0.1,
  max_time = 100,
  start_time = 0
)
```

## Arguments

- luck:

  Numeric between 0 and 1. Initial random quantile (luck).

- a_fun:

  Function of time .time returning the first distribution parameter
  (e.g., rate, shape, meanlog).

- b_fun:

  Function of time .time returning the second distribution parameter
  (e.g., scale, sdlog). Defaults to a function returning NA.

- dist:

  Character string specifying the distribution. Supported: "exp",
  "gamma", "lnorm", "norm", "weibull", "llogis", "gompertz".

- dt:

  Numeric. Time step increment to update parameters and survival.
  Default 0.1.

- max_time:

  Numeric. Max allowed event time to prevent infinite loops. Default
  100.

- start_time:

  Numeric. Time to use as a starting point of reference (e.g., curtime).

## Value

List with simulated time-to-event and final luck value.

## Details

The objective of this function is to avoid the user to have cycle events
with the only scope of updating some variables that depend on time and
re-evaluate a TTE. The idea is that this function should only be called
at start and when an event impacts a variable (e.g., stroke event
impacting death TTE), in which case it would need to be called again at
that point. In that case, the user would need to call e.g.,
`a <- qtimecov` with `max_time = curtime` arguments, and then call it
again with no max_time, and `luck = a$luck, start_time=a$tte` (so there
is no need to add curtime to the resulting time).

It's recommended to play with `dt` argument to balance running time and
precision of the estimates. For example, if we know we only update the
equation annually (not continuously), then we could just set `dt = 1`,
which would make computations faster.

## Examples

``` r
param_fun_factory <- function(p0, p1, p2, p3) {
  function(.time) p0 + p1*.time + p2*.time^2 + p3*(floor(.time) + 1)
}

set.seed(42)

# 1. Exponential Example
rate_exp <- param_fun_factory(0.1, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = rate_exp,
  dist = "exp"
)
#> $tte
#> [1] 24.62825
#> 
#> $luck
#> [1] 0.002820795
#> 


# 2. Gamma Example
shape_gamma <- param_fun_factory(2, 0, 0, 0)
rate_gamma <- param_fun_factory(0.2, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = shape_gamma,
  b_fun = rate_gamma,
  dist = "gamma"
)
#> $tte
#> [1] 22.32001
#> 
#> $luck
#> [1] 0.003261765
#> 


# 3. Lognormal Example
meanlog_lnorm <- param_fun_factory(log(10) - 0.5*0.5^2, 0, 0, 0)
sdlog_lnorm <- param_fun_factory(0.5, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = meanlog_lnorm,
  b_fun = sdlog_lnorm,
  dist = "lnorm"
)
#> $tte
#> [1] 6.655032
#> 
#> $luck
#> [1] 0.007686672
#> 


# 4. Normal Example
mean_norm <- param_fun_factory(10, 0, 0, 0)
sd_norm <- param_fun_factory(2, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = mean_norm,
  b_fun = sd_norm,
  dist = "norm"
)
#> $tte
#> [1] 11.9122
#> 
#> $luck
#> [1] 0.008791272
#> 


# 5. Weibull Example
shape_weibull <- param_fun_factory(2, 0, 0, 0)
scale_weibull <- param_fun_factory(10, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = shape_weibull,
  b_fun = scale_weibull,
  dist = "weibull"
)
#> $tte
#> [1] 10.13201
#> 
#> $luck
#> [1] 0.006391193
#> 


# 6. Loglogistic Example
shape_llogis <- param_fun_factory(2.5, 0, 0, 0)
scale_llogis <- param_fun_factory(7.6, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = shape_llogis,
  b_fun = scale_llogis,
  dist = "llogis"
)
#> $tte
#> [1] 7.836007
#> 
#> $luck
#> [1] 0.005926301
#> 


# 7. Gompertz Example
shape_gomp <- param_fun_factory(0.01, 0, 0, 0)
rate_gomp <- param_fun_factory(0.091, 0, 0, 0)
qtimecov(
  luck = runif(1),
  a_fun = shape_gomp,
  b_fun = rate_gomp,
  dist = "gompertz"
)
#> $tte
#> [1] 13.67996
#> 
#> $luck
#> [1] 0.008297281
#> 

#Time varying example, with change at time 8
rate_exp <- function(.time) 0.1 + 0.01*.time * 0.00001*.time^2
rate_exp2 <- function(.time) 0.2 + 0.02*.time
time_change <- 8
init_luck <- 0.95

a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
                      max_time = time_change)
qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
#> $tte
#> [1] 13.31393
#> 
#> $luck
#> [1] 0.001832564
#> 


#An example of how it would work in the model, this would also work with time varying covariates!
rate_exp <- function(.time) 0.1
rate_exp2 <- function(.time) 0.2
rate_exp3 <- function(.time) 0.3
time_change <- 10 #evt 1
time_change2 <- 15 #evt2
init_luck <- 0.95
#at start, we would just draw TTE
qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005)
#> $tte
#> [1] 29.95732
#> 
#> $luck
#> [1] 0.0002322466
#> 

#at event in which rate changes (at time 10) we need to do this:
a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
                      max_time = time_change)
new_luck <- a$luck
qtimecov(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
#> $tte
#> [1] 19.97866
#> 
#> $luck
#> [1] 0.0007320055
#> 

#at second  event in which rate changes again (at time 15) we need to do this:
a <- qtimecov(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005,
                      max_time = time_change2, start_time=a$tte)
new_luck <- a$luck
#final TTE is
qtimecov(luck = new_luck,a_fun = rate_exp3,dist = "exp", dt = 0.005, start_time=a$tte)
#> $tte
#> [1] 18.31911
#> 
#> $luck
#> [1] 0.001231515
#> 
```
