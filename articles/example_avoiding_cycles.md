# How to Avoid Using Cycles to Speed Up Model

## Introduction

This document shows how one can substantially reduce running time in
models where cycle events are used by utilizing three features of
WARDEN: 1) per cycle costs and qalys, 2) time to event predictions which
take into account time-varying covariates through
[`qtimecov()`](https://jsanchezalv.github.io/WARDEN/reference/qtimecov.md)
and 3) adjustment of ongoing outputs to be corrected by cycle-type of
adjustments (e.g., utility adjustment by age) through
[`adj_val()`](https://jsanchezalv.github.io/WARDEN/reference/adj_val.md).

Note as well that this approach to avoid cycles could smoothly handle
treatment waning, changes in covariates (like weight, BMI, etc), etc.

### Main options

``` r
library(WARDEN)
library(flexsurv)
#> Loading required package: survival
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(kableExtra)
#> 
#> Attaching package: 'kableExtra'
#> The following object is masked from 'package:dplyr':
#> 
#>     group_rows
library(purrr)
library(tidyr)
```

``` r
options(scipen = 999)
options(digits=3)
options(tibble.print_max = 50)
```

## Model using cycles

This is a very simple model, where we just have a time to death given by
an exponential, but the underlying risk changes every year given by a
patient-specific time covariate. This means that in the traditional
approach, we would create yearly cycles to make the adjustments to the
time to event.

Additionally, we have a multiplicative utility age adjustment factor, so
at every yearly cycle we also adjust for that factor.

Treatment costs are yearly as well, so we add them as an instantaneous
cost at every cycle, including time 0.

### Inputs

``` r

common_all_inputs <-add_item(input = {
  u_bs <- 0.8
  u_age <- 1 - seq(from = 0, to = 0.5, by = 0.005)
  bs_rate <- 0.07
  trt_eff <- 0.02
})  

common_pt_inputs <- add_item(input={
    bs_age <- rnorm(1,40,5)
    time_cov <- runif(1)/500
    luck <- runif(1)
}) 


unique_pt_inputs <- add_item(input = {
    cost_trt_ins <- ifelse(arm=="int",1000,500)
    bs_rate <- bs_rate - ifelse(arm=="int",trt_eff,0)
    bs_rate_f <- function(.time){bs_rate + (time_cov * floor(.time))}

})


init_event_list <- 
  add_tte(arm=c("int","noint"), evts = c("start","cycle","death") ,input={
    start <- 0
    cycle <- 1
    death <- qcond_exp(luck, rate = bs_rate)
  })



evt_react_list <-
  add_reactevt(name_evt = "start",
               input = {
                 new_rate <- bs_rate
                 q_default <- u_bs * u_age[min(100,floor(curtime + bs_age))]
                 cost_trt <- cost_trt_ins
               }) %>%
  add_reactevt(name_evt = "cycle",
               input = {
                 # if(curtime >= 100 - bs_age){
                 #   curtime <- Inf
                 # }
              
              new_event(c(cycle = curtime + 1))
                 
              q_default <- u_bs * u_age[min(100,floor(curtime + bs_age))]
              if(curtime<5){
              cost_trt <- cost_trt_ins
              }
              
              old_rate <- new_rate
              new_rate <- bs_rate_f(curtime)
              prev_surv <- 1 - pexp(curtime - 1, old_rate)
              cur_surv <- 1 - pexp(curtime, old_rate)
              
              luck <- luck_adj(prevsurv = prev_surv, cursurv = cur_surv, luck = luck, condq = TRUE)
              ttdeath <- qcond_exp(luck, rate = new_rate) + curtime
              modify_event(c(death = max(curtime, ttdeath)))

               })%>%
  add_reactevt(name_evt = "death",
               input = {
                 curtime <- Inf
               }) 

util_ongoing <- "q_default"

cost_instant <- "cost_trt"

results <- run_sim(  
  npats=5000,                               
  n_sim=1,                                  
  psa_bool = FALSE,                         
  arm_list = c("int", "noint"),             
  common_all_inputs = common_all_inputs,   
  common_pt_inputs = common_pt_inputs,      
  unique_pt_inputs = unique_pt_inputs,
  init_event_list = init_event_list,        
  evt_react_list = evt_react_list,          
  util_ongoing_list = util_ongoing,
  cost_instant_list = cost_instant,
  ipd = 2, input_out = c("luck","time_cov")
)
#> Analysis number: 1
#> Simulation number: 1
#> Patient-arm data aggregated across events by selecting the last value for input_out items.
#> Time to run simulation 1: 13.8s
#> Time to run analysis 1: 13.8s
#> Total time to run: 13.81s
#> Simulation finalized;
```

### Summary of Results

    #>                       int    noint
    #> costs             4298.77  2073.04
    #> dcosts               0.00  2225.73
    #> lys                 11.33     9.37
    #> dlys                 0.00     1.97
    #> qalys                6.87     5.74
    #> dqalys               0.00     1.13
    #> ICER                   NA  1131.35
    #> ICUR                   NA  1967.29
    #> INMB                   NA 54342.84
    #> costs_undisc      4543.20  2188.30
    #> dcosts_undisc        0.00  2354.90
    #> lys_undisc          16.26    12.60
    #> dlys_undisc          0.00     3.65
    #> qalys_undisc         9.61     7.59
    #> dqalys_undisc        0.00     2.02
    #> ICER_undisc            NA   644.34
    #> ICUR_undisc            NA  1165.14
    #> INMB_undisc            NA 98701.76
    #> cost_trt          4298.77  2073.04
    #> dcost_trt            0.00  2225.73
    #> cost_trt_undisc   4543.20  2188.30
    #> dcost_trt_undisc     0.00  2354.90
    #> luck                 0.03     0.04
    #> dluck                0.00    -0.01
    #> q_default            6.87     5.74
    #> dq_default           0.00     1.13
    #> q_default_undisc     9.61     7.59
    #> dq_default_undisc    0.00     2.02
    #> time_cov             0.00     0.00
    #> dtime_cov            0.00     0.00

| pat_id | arm   | total_lys | total_qalys | total_costs | total_lys_undisc | total_qalys_undisc | total_costs_undisc | cost_trt | q_default | q_default_undisc | cost_trt_undisc | nexttime | number_events |  luck | time_cov | simulation | sensitivity |
|-------:|:------|----------:|------------:|------------:|-----------------:|-------------------:|-------------------:|---------:|----------:|-----------------:|----------------:|---------:|--------------:|------:|---------:|-----------:|------------:|
|      1 | int   |      8.73 |        5.44 |        4717 |            10.10 |               6.28 |               5000 |     4717 |      5.44 |             6.28 |            5000 |     75.2 |            12 | 0.005 |    0.000 |          1 |           1 |
|      1 | noint |      6.56 |        4.12 |        2359 |             7.29 |               4.57 |               2500 |     2359 |      4.12 |             4.57 |            2500 |     42.6 |             9 | 0.020 |    0.000 |          1 |           1 |
|      2 | int   |     17.85 |       11.09 |        4717 |            25.37 |              15.61 |               5000 |     4717 |     11.09 |            15.61 |            5000 |    375.7 |            27 | 0.020 |    0.000 |          1 |           1 |
|      2 | noint |     14.28 |        9.03 |        2359 |            18.55 |              11.66 |               2500 |     2359 |      9.03 |            11.66 |            2500 |    208.1 |            20 | 0.040 |    0.000 |          1 |           1 |
|      3 | int   |      4.74 |        2.98 |        4717 |             5.11 |               3.21 |               5000 |     4717 |      2.98 |             3.21 |            5000 |     25.2 |             7 | 0.006 |    0.000 |          1 |           1 |
|      3 | noint |      3.47 |        2.19 |        1914 |             3.66 |               2.31 |               2000 |     1914 |      2.19 |             2.31 |            2000 |     13.3 |             5 | 0.045 |    0.000 |          1 |           1 |
|      4 | int   |     12.88 |        8.09 |        4717 |            16.21 |              10.14 |               5000 |     4717 |      8.09 |            10.14 |            5000 |    168.4 |            18 | 0.011 |    0.000 |          1 |           1 |
|      4 | noint |      9.85 |        6.27 |        2359 |            11.64 |               7.39 |               2500 |     2359 |      6.27 |             7.39 |            2500 |     89.3 |            13 | 0.045 |    0.000 |          1 |           1 |
|      5 | int   |     20.75 |       11.95 |        4717 |            32.15 |              18.19 |               5000 |     4717 |     11.95 |            18.19 |            5000 |    592.3 |            34 | 0.012 |    0.001 |          1 |           1 |
|      5 | noint |     17.94 |       10.50 |        2359 |            25.56 |              14.79 |               2500 |     2359 |     10.50 |            14.79 |            2500 |    376.1 |            27 | 0.051 |    0.001 |          1 |           1 |

## Model avoiding cycles

We now avoid to use cycles.

To fix the time to death with time varying covariates, we use `timecov`
with yearly updates (`dt = 1`) and passing the parameter as a function
of time, so the function automatically iterates over each year to obtain
the relevant time to event. This means we can obtain the final time to
event in one go at start. While still iterative, it’s less time
consuming than having events created at every year. Note that if we had
a continuous time changing covariate, we could still approximate it
quite well by using `dt = 0.01` (every 3.5 days updates).

The multiplicative utility age adjustment factor is implemented by
[`adj_val()`](https://jsanchezalv.github.io/WARDEN/reference/adj_val.md),
a function that given an expression that changes with time, reweighs the
evaluated values by the ongoing discounting factor that would apply so
that even if a single value is returned, the discounted outcome will
still be the same as in the traditional cycle approach. Note however
that undiscounted outcomes will not be correct (in this case we do not
care about it, but we could just create another value like
`other_q_default` where the
[`adj_val()`](https://jsanchezalv.github.io/WARDEN/reference/adj_val.md)
is used with 0 discounting, and that value will be correct
undiscounted). The very small difference INMB comes from rounding
effects.

Treatment costs are yearly, so we just use the per cycle approach
instead of instantaneous costs at every year.

### Inputs

``` r

common_all_inputs <-add_item(input = {
  u_bs <- 0.8
  u_age <- 1 - seq(from = 0, to = 0.5, by = 0.005)
  bs_rate <- 0.07
  trt_eff <- 0.02
  
  #set costs as cycle
  cost_trt_cycle_l <- 1
  cost_trt_cycle_starttime <- 0
  cost_trt_max_cycles <- 5
})  

common_pt_inputs <- add_item(input={
    bs_age <- rnorm(1,40,5)
    time_cov <- runif(1)/500
    luck <- runif(1)
}) 


unique_pt_inputs <- add_item(input = {
    cost_trt <- ifelse(arm=="int",1000,500)
    bs_rate <- bs_rate - ifelse(arm=="int",trt_eff,0)
    #rate is a function of time, with yearly changes
    bs_rate_f <- function(.time){bs_rate + (time_cov * floor(.time))}
   
})


init_event_list <- 
  add_tte(arm=c("int","noint"), evts = c("start","death") ,input={
    start <- 0
    #we can immediately know what the final TTE will be even with time changing covariates
    death <- qtimecov(luck = luck,a_fun = bs_rate_f,dist = "exp", dt = 1)$tte
  })



evt_react_list <-
  add_reactevt(name_evt = "start",
               input = {
                 new_rate <- bs_rate
                 #this is the way to obtain utility adjustment by age as a single value, it reweighs the u_age by their discounting value so final discounted outcomes are correct
                 adj_factor <- adj_val(
                   curtime,
                   next_event()$time,
                   by = 1,
                   u_age[min(100,floor(.time + bs_age))],
                   discount = drq
                   )
                 q_default <- u_bs * adj_factor
               }) %>%
  add_reactevt(name_evt = "death",
               input = {
                 curtime <- Inf
               }) 

util_ongoing <- "q_default"

cost_cycle <- "cost_trt"

results2 <- run_sim(  
  npats=5000,                               
  n_sim=1,                                  
  psa_bool = FALSE,                         
  arm_list = c("int", "noint"),             
  common_all_inputs = common_all_inputs,   
  common_pt_inputs = common_pt_inputs,      
  unique_pt_inputs = unique_pt_inputs,
  init_event_list = init_event_list,        
  evt_react_list = evt_react_list,          
  util_ongoing_list = util_ongoing,
  cost_cycle_list = cost_cycle,
  ipd = 2, input_out = c("luck","time_cov")
)
#> Analysis number: 1
#> Simulation number: 1
#> Patient-arm data aggregated across events by selecting the last value for input_out items.
#> Time to run simulation 1: 3.92s
#> Time to run analysis 1: 3.92s
#> Total time to run: 3.92s
#> Simulation finalized;
```

### Summary of Results

    #>                               int     noint
    #> costs                     4298.77   2073.04
    #> dcosts                       0.00   2225.73
    #> lys                         11.33      9.37
    #> dlys                         0.00      1.97
    #> qalys                        6.87      5.74
    #> dqalys                       0.00      1.13
    #> ICER                           NA   1131.41
    #> ICUR                           NA   1967.37
    #> INMB                           NA  54340.53
    #> costs_undisc              4543.20   2188.30
    #> dcosts_undisc                0.00   2354.90
    #> lys_undisc                  16.26     12.60
    #> dlys_undisc                  0.00      3.65
    #> qalys_undisc                 9.78      7.68
    #> dqalys_undisc                0.00      2.10
    #> ICER_undisc                    NA    644.83
    #> ICUR_undisc                    NA   1120.61
    #> INMB_undisc                    NA 102717.14
    #> cost_trt                  4298.77   2073.04
    #> dcost_trt                    0.00   2225.73
    #> cost_trt_cycle_l             2.00      2.00
    #> dcost_trt_cycle_l            0.00      0.00
    #> cost_trt_cycle_starttime     0.00      0.00
    #> dcost_trt_cycle_starttime    0.00      0.00
    #> cost_trt_max_cycles         10.00     10.00
    #> dcost_trt_max_cycles         0.00      0.00
    #> cost_trt_undisc           4543.20   2188.30
    #> dcost_trt_undisc             0.00   2354.90
    #> luck                         0.50      0.50
    #> dluck                        0.00      0.00
    #> q_default                    6.87      5.74
    #> dq_default                   0.00      1.13
    #> q_default_undisc             9.78      7.68
    #> dq_default_undisc            0.00      2.10
    #> time_cov                     0.00      0.00
    #> dtime_cov                    0.00      0.00

| pat_id | arm   | total_lys | total_qalys | total_costs | total_lys_undisc | total_qalys_undisc | total_costs_undisc | cost_trt | cost_trt_cycle_l | cost_trt_cycle_starttime | cost_trt_max_cycles | q_default | q_default_undisc | cost_trt_undisc | nexttime | number_events |  luck | time_cov | simulation | sensitivity |
|-------:|:------|----------:|------------:|------------:|-----------------:|-------------------:|-------------------:|---------:|-----------------:|-------------------------:|--------------------:|----------:|-----------------:|----------------:|---------:|--------------:|------:|---------:|-----------:|------------:|
|      1 | int   |      8.73 |        5.44 |        4717 |            10.10 |               6.29 |               5000 |     4717 |                2 |                        0 |                  10 |      5.44 |             6.29 |            5000 |    20.20 |             2 | 0.403 |    0.000 |          1 |           1 |
|      1 | noint |      6.56 |        4.12 |        2359 |             7.29 |               4.58 |               2500 |     2359 |                2 |                        0 |                  10 |      4.12 |             4.58 |            2500 |    14.57 |             2 | 0.403 |    0.000 |          1 |           1 |
|      2 | int   |     17.85 |       11.09 |        4717 |            25.37 |              15.77 |               5000 |     4717 |                2 |                        0 |                  10 |     11.09 |            15.77 |            5000 |    50.74 |             2 | 0.736 |    0.000 |          1 |           1 |
|      2 | noint |     14.28 |        9.03 |        2359 |            18.55 |              11.73 |               2500 |     2359 |                2 |                        0 |                  10 |      9.03 |            11.73 |            2500 |    37.10 |             2 | 0.736 |    0.000 |          1 |           1 |
|      3 | int   |      4.74 |        2.98 |        4717 |             5.11 |               3.21 |               5000 |     4717 |                2 |                        0 |                  10 |      2.98 |             3.21 |            5000 |    10.22 |             2 | 0.226 |    0.000 |          1 |           1 |
|      3 | noint |      3.47 |        2.19 |        1914 |             3.66 |               2.31 |               2000 |     1914 |                2 |                        0 |                  10 |      2.19 |             2.31 |            2000 |     7.32 |             2 | 0.226 |    0.000 |          1 |           1 |
|      4 | int   |     12.88 |        8.09 |        4717 |            16.21 |              10.18 |               5000 |     4717 |                2 |                        0 |                  10 |      8.09 |            10.18 |            5000 |    32.41 |             2 | 0.560 |    0.000 |          1 |           1 |
|      4 | noint |      9.85 |        6.27 |        2359 |            11.64 |               7.41 |               2500 |     2359 |                2 |                        0 |                  10 |      6.27 |             7.41 |            2500 |    23.29 |             2 | 0.560 |    0.000 |          1 |           1 |
|      5 | int   |     20.75 |       11.95 |        4717 |            32.15 |              18.51 |               5000 |     4717 |                2 |                        0 |                  10 |     11.95 |            18.51 |            5000 |    64.31 |             2 | 0.877 |    0.001 |          1 |           1 |
|      5 | noint |     17.94 |       10.50 |        2359 |            25.56 |              14.96 |               2500 |     2359 |                2 |                        0 |                  10 |     10.50 |            14.96 |            2500 |    51.12 |             2 | 0.877 |    0.001 |          1 |           1 |
