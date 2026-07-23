# How to Use the Automatic Input Selector

``` r

library("dplyr")
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library("kableExtra")
#> 
#> Attaching package: 'kableExtra'
#> The following object is masked from 'package:dplyr':
#> 
#>     group_rows
library("knitr")
library("purrr")
library("MASS")
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
library("WARDEN")
```

``` r

options(scipen = 999)
options(digits=3)
options(tibble.print_max = 50)
```

## Introduction

This document explains how to use the set of functions related to
automatic input selector, particularly
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md),
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md),
[`pick_psa()`](https://jsanchezalv.github.io/WARDEN/reference/pick_psa.md),
and
[`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md),
so that the model takes care of everything in terms of running a
deterministic analysis, DSA, PSA, probabilistic DSA, or scenario
analyses.

We use
[`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md)
for all the examples below, in previous versions of WARDEN we
distinguished between add_item and add_item2, but now the behavior has
been integrated into a single function that can handle both approaches
of adding inputs.

### In a Nutshell

The key function
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md).
This function essentially hides a loop, which iterates over each of the
inputs and depending on whether the parameter is a vector or single
length, and whether the PSA or scenario flags are active, proceeds to
select the right value. It also can consider multiple parameters being
covaried at the same time in a scenario or DSA analysis.

The
[`pick_psa()`](https://jsanchezalv.github.io/WARDEN/reference/pick_psa.md)
function is a wrapper that just calls the corresponding function from
its first argument, so it can be used to draw from distributions. The
recommendation is to use the “r” functions like
[`rnorm()`](https://rdrr.io/r/stats/Normal.html),
[`runif()`](https://rdrr.io/r/stats/Uniform.html), etc. The random seed
is handled automatically by the model to ensure that inputs are drawn
appropiately.

Finally, the
[`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md)
is a function that will generate a vector of 0 and 1 and is used for
scenario analyses and DSA, as e.g., in a DSA we need to iterate over the
corresponding parameters.

The function
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md)
is essentially a wrapper for
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
which will make sure to call the correct parameters, sensitivity
analysis required, etc.

For now, we can start with a simple example, see the data below, where
we have a set of parameters, with some base case values, PSA parameters,
DSA, scenario, and whether the parameter would be active on a PSA
(vector of 0 and 1s):

| parameter_name | base_value | PSA_dist | a | b | n | DSA_min | DSA_max | scenario_1 | scenario_2 | psa_indicators |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| util.sick | 0.8 | rnorm | 0.8 | 0.16 | 1 | 0.6 | 0.9 | 0.6 | 0.9 | 1 |
| util.sicker | 0.5 | rbeta_mse | 0.5 | 0.1 | 1 | 0.3 | 0.7 | 0.3 | 0.7 | 1 |
| cost.sick | 3000 | rgamma_mse | 3000 | 600 | 1 | 1000 | 5000 | 1000 | 5000 | 1 |
| cost.sicker | 7000 | rgamma_mse | 7000 | 1400 | 1 | 5000 | 9000 | 5000 | 9000 | 1 |
| cost.int | 1000 | rgamma_mse | 1000 | 200 | 1 | 800 | 2000 | 800 | 2000 | 0 |
| coef_noint | -1.6094379124341 | rnorm | -1.6094379124341 | 0.32188758248682 | 1 | -2.30258509299405 | -0.916290731874155 | -2.30258509299405 | -0.916290731874155 | 0 |
| HR_int | 0.8 | rlnorm | -0.22314355131421 | 0.0446287102628419 | 1 | 0.5 | 0.9 | 0.5 | 0.9 | 0 |

### The basic setup

We’ll build our case step by step. in
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md),
we need to set which values would be the base case (`base`), the PSA
(`psa`), and the sensitivity values (`sens`). We also need to set
whether we are currently in a PSA or not (`psa_ind`), and whether we are
in a sensitivity analysis (`sens_ind`). Finally, we need to give some
guidance if we are in a sensitivity analysis on which parameters to vary
at this iteration (i.e., if it’s iterating over the minimum range of the
DSA, is it `util.sick`? is it `util.sicker`?) and to give the names of
the parameters to export the named list. The simplest case is when we
don’t have DSA or scenario analysis. In that case the function is
relatively straightforward.

We set our base case values, we set the PSA values using the parameters
in our data, indicate that we are indeed in a PSA analysis and NOT in a
scenario analysis (we could have both of them being `TRUE`), we indicate
that if it was a scenario analysis all parameters would be included, and
finally we indicate that for the PSA, the first 4 will be included and
the remaining 3 will not (defaulting to the base values).

``` r


  pick_val_v(base        = l_inputs[["base_value"]],
             psa         = pick_psa(
                 l_inputs[["PSA_dist"]],
                 l_inputs[["n"]],
                 l_inputs[["a"]],
                 l_inputs[["b"]]),
             psa_ind     = TRUE, #PSA is active
             sens_ind    = FALSE, #No scenario analysis
             names_out   = l_inputs[["parameter_name"]],
             indicator   = rep(1,7), #This is only relevant for scenario analysis, set to all 1s
             indicator_psa = l_inputs[["psa_indicators"]],#vector of 4 1s and 3 0s.
             deploy_env = FALSE 
  ) 
#> $util.sick
#> [1] 0.576
#> 
#> $util.sicker
#> [1] 0.529
#> 
#> $cost.sick
#> [1] 1671
#> 
#> $cost.sicker
#> [1] 6114
#> 
#> $cost.int
#> [1] 1000
#> 
#> $coef_noint
#> [1] -1.61
#> 
#> $HR_int
#> [1] 0.8
```

See below a more complex example, now including scenario analyses/DSA,
and letting the model handling more things automatically.

Note that
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
has a few arguments set to variables: `sens_name_used`, `psa_bool`,
`sensitivity_bool` and `indicators`. All except `indicators` are set by
the model automatically, and are obtained from
[`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
(from `sensitivity_names`, `psa_bool` and `sensitivity_bool`,
respectively). `sens_name_used` is selected automatically by the model,
so for example if `sensitivity_names = "DSA_min", "DSA_max"`, the model
automatically knows which one is currently active, so it will start with
`"DSA_min"`, and once it has iterated through all the relevant
parameters, it will go to `"DSA_max"`. `psa_bool` and `sensitivity_bool`
are boolean flags for whether the current analysis is
deterministic/probabilistic and standard/scenario(or DSA). Since in this
example we are not using `run_sim`, we are setting this in advance.

`indicators` is an object that we need to set when we are setting our
inputs with
[`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md),
and it will be used only for scenario/DSA analyses (so when
`send_ind = sensitivity_bool = TRUE`). In this case, we set it to all 1s
(with length equal to the number of vectors) whenever we are not in a
scenario analysis, and otherwise we use the `create_indicators` function
that takes a few extra arguments. Alternatively, the user does not need
to worry about the indicators as it will be taken care by the
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md)
function (see below).

If the user wants to use the legacy option with
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
instead of
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md),
then
[`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md)
creates a vector of 0 and 1s, taking value 1 at the right index that is
going to be varied. It takes a few extra objects that are given
automatically by the model: `sens`, `n_sensitivity`, and
`sensitivity_names`. `sensitivity_names` comes from `run_sim` as
mentioned above, and `n_sensitivity` also must be declared in `run_sim`,
and is the number of parameters to be varied, in this case, 7). `sens`
is the index of the analysis currently being run, in this case we assume
we start with the first index so it will return a vector of 1 followed
by 6 0s.

Note that we do not necessarily have to use the `indicator_psa` argument
of
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md).
If left as is, it will understand that all parameters would be included
in the PSA. If not, we would need to provide that argument declaring
which parameters are excluded or included (vector of 0 and 1s).

To recap: in this case, we artificially set some key variables
beforehand, and then we execute the
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
function. We first run it as a deterministic analysis.

``` r


sens <- 1
n_sensitivity <- length(l_inputs[[1]])
sensitivity_names <- c("DSA_min", "DSA_max")
#Vector of length 7, a 1 followed by 6 0s, if sens was 2 it would be a c(0,1,0,0,0,0,0) vector 
create_indicators(sens,n_sensitivity*length(sensitivity_names),rep(1,length(l_inputs[[1]])))
#> [1] 1 0 0 0 0 0 0

sens_name_used <- "DSA_min"
psa_bool <- FALSE
sensitivity_bool <- FALSE

#deterministic
indicators <-  if(sensitivity_bool){ create_indicators(sens,n_sensitivity*length(sensitivity_names),rep(1,length(l_inputs[[1]])))}else{
                                rep(1,length(l_inputs[[1]]))} 

#DETERMINISTIC
as.data.frame(
  pick_val_v(base        = l_inputs[["base_value"]],
           psa         = pick_psa(
             l_inputs[["PSA_dist"]],
             l_inputs[["n"]],
             l_inputs[["a"]],
             l_inputs[["b"]]),
           sens        = l_inputs[[sens_name_used]], #e.g., sens_name_used = "DSA_min"
           psa_ind     = psa_bool, #FALSE
           sens_ind    = sensitivity_bool, #FALSE
           indicator   = indicators, #all 1s, or a vector of 1 1 and the rest 0s.
           names_out   = l_inputs[["parameter_name"]],
           indicator_psa = l_inputs[["psa_indicators"]],
           deploy_env = FALSE
           )
) 
#>   util.sick util.sicker cost.sick cost.sicker cost.int coef_noint HR_int
#> 1       0.8         0.5      3000        7000     1000      -1.61    0.8
```

Now it’s very easy to switch to probabilistic analysis.

``` r

#PSA
psa_bool <- TRUE
as.data.frame(
  pick_val_v(base        = l_inputs[["base_value"]],
           psa         = pick_psa(
             l_inputs[["PSA_dist"]],
             l_inputs[["n"]],
             l_inputs[["a"]],
             l_inputs[["b"]]),
           sens        = l_inputs[[sens_name_used]], #e.g., sens_name_used = "DSA_min"
           psa_ind     = psa_bool, #FALSE
           sens_ind    = sensitivity_bool, #FALSE
           indicator   = indicators, #all 1s, or a vector of 1 1 and the rest 0s.
           names_out   = l_inputs[["parameter_name"]],
           indicator_psa = l_inputs[["psa_indicators"]] ,
           deploy_env = FALSE
           )
)
#>   util.sick util.sicker cost.sick cost.sicker cost.int coef_noint HR_int
#> 1     0.761       0.467      2620        5725     1000      -1.61    0.8
```

And it’s very easy to switch to probabilistic scenario analysis.

``` r

#Probabilistic DSA, first parameter being varied as sens = 1
psa_bool <- TRUE
sensitivity_bool <- TRUE
indicators <-  if(sensitivity_bool){ create_indicators(sens,n_sensitivity*length(sensitivity_names),rep(1,length(l_inputs[[1]])))}else{
                                rep(1,length(l_inputs[[1]]))} 
as.data.frame(
  pick_val_v(base        = l_inputs[["base_value"]],
           psa         = pick_psa(
             l_inputs[["PSA_dist"]],
             l_inputs[["n"]],
             l_inputs[["a"]],
             l_inputs[["b"]]),
           sens        = l_inputs[[sens_name_used]], #e.g., sens_name_used = "DSA_min"
           psa_ind     = psa_bool, #FALSE
           sens_ind    = sensitivity_bool, #FALSE
           indicator   = indicators, #all 1s, or a vector of 1 1 and the rest 0s.
           names_out   = l_inputs[["parameter_name"]],
           indicator_psa = l_inputs[["psa_indicators"]] ,
           deploy_env = FALSE
           )
) 
#>   util.sick util.sicker cost.sick cost.sicker cost.int coef_noint HR_int
#> 1       0.6        0.48      3989        5818     1000      -1.61    0.8
```

### Integrating into a model

We now use this setup in a very simple model, where we run a
deterministic, probabilistic, probabilistic DSA and deterministic
scenario analysis. With
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md),
the `n_sensitivity` and `sensitivity_names` for DSA are auto-detected
from the block’s metadata — no manual `i_sensitivity` or `sens_iterator`
setup is needed. For scenario analyses, `sensitivity_names` must still
be provided explicitly in
[`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
to understand which sensitivity analysis should be run, but
`n_sensitivity` is auto-calculated.

``` r


rm(sens, sens_name_used, sensitivity_bool, psa_bool) #remove global objects that may confuse program

i_simple <- input_block(
  base           = l_inputs[["base_value"]],
  psa            = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
  sens           = l_inputs,
  names_out      = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]],
  dsa_names      = c("DSA_min", "DSA_max"))

i_arm <- add_item(q_default = util.sick,
         c_default = cost.sick + if(arm=="int"){cost.int}else{0})


init_event_list <- 
  add_tte(arm=c("noint","int"), evts = c("a1","b1") ,input={
    a1 <- 0
    b1 <- 2
  })

evt_react_list <-
  add_reactevt(name_evt = "a1",
               input = {})  |>
  add_reactevt(name_evt = "b1",
               input = {
                 q_default = 0
                 c_default = 0 
                 curtime = Inf
               }) 

util_ongoing <- "q_default"
cost_ongoing <- "c_default"

results <- run_sim(  
  npats=5,                               
  n_sim=1,                                  
  psa_bool = FALSE,                        
  arm_list = c("int", "noint"),             
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  init_event_list = init_event_list,        
  evt_react_list = evt_react_list,          
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1
)
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.05s
#> Time to run analysis 1: 0.05s
#> Total time to run: 0.06s
#> Simulation finalized;

summary_results_sim(results[[1]])  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

|                   | int                  | noint                   |
|:------------------|:---------------------|:------------------------|
| costs             | 7,768 (7,768; 7,768) | 5,826 (5,826; 5,826)    |
| dcosts            | 0 (0; 0)             | 1,942 (1,942; 1,942)    |
| lys               | 1.94 (1.94; 1.94)    | 1.94 (1.94; 1.94)       |
| dlys              | 0 (0; 0)             | 0 (0; 0)                |
| qalys             | 1.55 (1.55; 1.55)    | 1.55 (1.55; 1.55)       |
| dqalys            | 0 (0; 0)             | 0 (0; 0)                |
| ICER              | NaN (NA; NA)         | Inf (Inf; Inf)          |
| ICUR              | NaN (NA; NA)         | Inf (Inf; Inf)          |
| INMB              | NaN (NA; NA)         | -1,942 (-1,942; -1,942) |
| costs_undisc      | 8,000 (8,000; 8,000) | 6,000 (6,000; 6,000)    |
| dcosts_undisc     | 0 (0; 0)             | 2,000 (2,000; 2,000)    |
| lys_undisc        | 2 (2; 2)             | 2 (2; 2)                |
| dlys_undisc       | 0 (0; 0)             | 0 (0; 0)                |
| qalys_undisc      | 1.6 (1.6; 1.6)       | 1.6 (1.6; 1.6)          |
| dqalys_undisc     | 0 (0; 0)             | 0 (0; 0)                |
| ICER_undisc       | NaN (NA; NA)         | Inf (Inf; Inf)          |
| ICUR_undisc       | NaN (NA; NA)         | Inf (Inf; Inf)          |
| INMB_undisc       | NaN (NA; NA)         | -2,000 (-2,000; -2,000) |
| c_default         | 7,768 (7,768; 7,768) | 5,826 (5,826; 5,826)    |
| dc_default        | 0 (0; 0)             | 1,942 (1,942; 1,942)    |
| c_default_undisc  | 8,000 (8,000; 8,000) | 6,000 (6,000; 6,000)    |
| dc_default_undisc | 0 (0; 0)             | 2,000 (2,000; 2,000)    |
| q_default         | 1.55 (1.55; 1.55)    | 1.55 (1.55; 1.55)       |
| dq_default        | 0 (0; 0)             | 0 (0; 0)                |
| q_default_undisc  | 1.6 (1.6; 1.6)       | 1.6 (1.6; 1.6)          |
| dq_default_undisc | 0 (0; 0)             | 0 (0; 0)                |

``` r


results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1
)
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.05s
#> Simulation number: 2
#> Time to run simulation 2: 0.05s
#> Time to run analysis 1: 0.1s
#> Total time to run: 0.11s
#> Simulation finalized;

summary_results_sim(results[[1]])  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

|                   | int                    | noint                   |
|:------------------|:-----------------------|:------------------------|
| costs             | 9,783 (7,895; 11,671)  | 7,841 (5,953; 9,729)    |
| dcosts            | 0 (0; 0)               | 1,942 (1,942; 1,942)    |
| lys               | 1.94 (1.94; 1.94)      | 1.94 (1.94; 1.94)       |
| dlys              | 0 (0; 0)               | 0 (0; 0)                |
| qalys             | 1.16 (1.02; 1.29)      | 1.16 (1.02; 1.29)       |
| dqalys            | 0 (0; 0)               | 0 (0; 0)                |
| ICER              | NaN (NA; NA)           | Inf (Inf; Inf)          |
| ICUR              | NaN (NA; NA)           | Inf (Inf; Inf)          |
| INMB              | NaN (NA; NA)           | -1,942 (-1,942; -1,942) |
| costs_undisc      | 10,075 (8,131; 12,020) | 8,075 (6,131; 10,020)   |
| dcosts_undisc     | 0 (0; 0)               | 2,000 (2,000; 2,000)    |
| lys_undisc        | 2 (2; 2)               | 2 (2; 2)                |
| dlys_undisc       | 0 (0; 0)               | 0 (0; 0)                |
| qalys_undisc      | 1.19 (1.06; 1.33)      | 1.19 (1.06; 1.33)       |
| dqalys_undisc     | 0 (0; 0)               | 0 (0; 0)                |
| ICER_undisc       | NaN (NA; NA)           | Inf (Inf; Inf)          |
| ICUR_undisc       | NaN (NA; NA)           | Inf (Inf; Inf)          |
| INMB_undisc       | NaN (NA; NA)           | -2,000 (-2,000; -2,000) |
| c_default         | 9,783 (7,895; 11,671)  | 7,841 (5,953; 9,729)    |
| dc_default        | 0 (0; 0)               | 1,942 (1,942; 1,942)    |
| c_default_undisc  | 10,075 (8,131; 12,020) | 8,075 (6,131; 10,020)   |
| dc_default_undisc | 0 (0; 0)               | 2,000 (2,000; 2,000)    |
| q_default         | 1.16 (1.02; 1.29)      | 1.16 (1.02; 1.29)       |
| dq_default        | 0 (0; 0)               | 0 (0; 0)                |
| q_default_undisc  | 1.19 (1.06; 1.33)      | 1.19 (1.06; 1.33)       |
| dq_default_undisc | 0 (0; 0)               | 0 (0; 0)                |

``` r

# DSA analyses — n_sensitivity and sensitivity_names auto-detected from input_block() metadata
results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1,
  sensitivity_bool = TRUE,
  input_out = unlist(l_inputs[["parameter_name"]])
)
#> n_sensitivity auto-detected as 7 from input_block metadata
#> sensitivity_names auto-set to c("DSA_min", "DSA_max") from input_block metadata
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.08s
#> Time to run analysis 1: 0.15s
#> Analysis number: 2
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 2: 0.12s
#> Analysis number: 3
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 3: 0.12s
#> Analysis number: 4
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.08s
#> Time to run analysis 4: 0.14s
#> Analysis number: 5
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 5: 0.15s
#> Analysis number: 6
#> Simulation number: 1
#> Time to run simulation 1: 0.08s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 6: 0.15s
#> Analysis number: 7
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 7: 0.13s
#> Analysis number: 8
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 8: 0.12s
#> Analysis number: 9
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 9: 0.12s
#> Analysis number: 10
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 10: 0.12s
#> Analysis number: 11
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 11: 0.12s
#> Analysis number: 12
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 12: 0.12s
#> Analysis number: 13
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 13: 0.12s
#> Analysis number: 14
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 14: 0.13s
#> Total time to run: 1.81s
#> Simulation finalized;

summary_results_sens(results)
#>          arm analysis analysis_name     variable                 value
#>       <char>    <int>        <char>       <fctr>                <char>
#>    1:    int        1       DSA_min        costs 9,783 (7,895; 11,671)
#>    2:  noint        1       DSA_min        costs  7,841 (5,953; 9,729)
#>    3:    int        2       DSA_min        costs 9,783 (7,895; 11,671)
#>    4:  noint        2       DSA_min        costs  7,841 (5,953; 9,729)
#>    5:    int        3       DSA_min        costs  3,884 (3,884; 3,884)
#>   ---                                                                 
#> 1116:  noint       12       DSA_max dutil.sicker              0 (0; 0)
#> 1117:    int       13       DSA_max dutil.sicker              0 (0; 0)
#> 1118:  noint       13       DSA_max dutil.sicker              0 (0; 0)
#> 1119:    int       14       DSA_max dutil.sicker              0 (0; 0)
#> 1120:  noint       14       DSA_max dutil.sicker              0 (0; 0)

data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))
data_sensitivity |> group_by(sensitivity) |> summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),mean)  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| sensitivity | util.sick | util.sicker | cost.sick | cost.sicker | cost.int | coef_noint | HR_int |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.600 | 0.554 | 4038 | 8633 | 1000 | -1.609 | 0.8 |
| 2 | 0.597 | 0.300 | 4038 | 8633 | 1000 | -1.609 | 0.8 |
| 3 | 0.597 | 0.554 | 1000 | 8633 | 1000 | -1.609 | 0.8 |
| 4 | 0.597 | 0.554 | 4038 | 5000 | 1000 | -1.609 | 0.8 |
| 5 | 0.597 | 0.554 | 4038 | 8633 | 800 | -1.609 | 0.8 |
| 6 | 0.597 | 0.554 | 4038 | 8633 | 1000 | -2.303 | 0.8 |
| 7 | 0.597 | 0.554 | 4038 | 8633 | 1000 | -1.609 | 0.5 |
| 8 | 0.900 | 0.554 | 4038 | 8633 | 1000 | -1.609 | 0.8 |
| 9 | 0.597 | 0.700 | 4038 | 8633 | 1000 | -1.609 | 0.8 |
| 10 | 0.597 | 0.554 | 5000 | 8633 | 1000 | -1.609 | 0.8 |
| 11 | 0.597 | 0.554 | 4038 | 9000 | 1000 | -1.609 | 0.8 |
| 12 | 0.597 | 0.554 | 4038 | 8633 | 2000 | -1.609 | 0.8 |
| 13 | 0.597 | 0.554 | 4038 | 8633 | 1000 | -0.916 | 0.8 |
| 14 | 0.597 | 0.554 | 4038 | 8633 | 1000 | -1.609 | 0.9 |

``` r

# Scenario analyses — n_sensitivity auto-detected as 1; sensitivity_names must be provided
results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1,
  sensitivity_bool = TRUE,
  sensitivity_names = c("scenario_1", "scenario_2"),
  input_out = unlist(l_inputs[["parameter_name"]])
)
#> n_sensitivity auto-detected as 7 from input_block metadata
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 1: 0.12s
#> Analysis number: 2
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 2: 0.12s
#> Total time to run: 0.24s
#> Simulation finalized;

summary_results_sens(results)
#>         arm analysis analysis_name     variable                   value
#>      <char>    <int>        <char>       <fctr>                  <char>
#>   1:    int        1    scenario_1        costs    3,496 (3,496; 3,496)
#>   2:  noint        1    scenario_1        costs    1,942 (1,942; 1,942)
#>   3:    int        2    scenario_2        costs 13,594 (13,594; 13,594)
#>   4:  noint        2    scenario_2        costs    9,710 (9,710; 9,710)
#>   5:    int        1    scenario_1       dcosts                0 (0; 0)
#>  ---                                                                   
#> 156:  noint        2    scenario_2  util.sicker          0.7 (0.7; 0.7)
#> 157:    int        1    scenario_1 dutil.sicker                0 (0; 0)
#> 158:  noint        1    scenario_1 dutil.sicker                0 (0; 0)
#> 159:    int        2    scenario_2 dutil.sicker                0 (0; 0)
#> 160:  noint        2    scenario_2 dutil.sicker                0 (0; 0)

data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))


data_sensitivity |> group_by(sensitivity) |> summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),mean)  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| sensitivity | util.sick | util.sicker | cost.sick | cost.sicker | cost.int | coef_noint | HR_int |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.6 | 0.3 | 1000 | 5000 | 800 | -2.303 | 0.5 |
| 2 | 0.9 | 0.7 | 5000 | 9000 | 2000 | -0.916 | 0.9 |

**Alternative: legacy approach using `pick_val_v()`**

Here we provide the legacy approach to the same model run.

``` r

# Build sensitivity indicator object and pick_val_v-based input block
i_sensitivity <- add_item(
  iterator_sensitivity = sens_iterator(sens, n_sensitivity)
) |> add_item(
  indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min", "DSA_max")) {
    create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                      rep(1, length(l_inputs[[1]])))
  } else { rep(1, length(l_inputs[[1]])) }
)

i_simple <- add_item() |>
  add_item(
    pick_val_v(
      base          = l_inputs[["base_value"]],
      psa           = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
      sens          = l_inputs[[sens_name_used]],
      psa_ind       = psa_bool,
      sens_ind      = sensitivity_bool,
      indicator     = indicators,
      names_out     = l_inputs[["parameter_name"]],
      indicator_psa = l_inputs[["psa_indicators"]]
    )
  )

# PSA run
run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
        common_all_inputs=i_simple, unique_pt_inputs=i_arm,
        init_event_list=init_event_list, evt_react_list=evt_react_list,
        util_ongoing_list=util_ongoing, cost_ongoing_list=cost_ongoing, ipd=1,
        sensitivity_inputs=i_sensitivity, sensitivity_names=NULL,
        sensitivity_bool=FALSE, n_sensitivity=1)

# DSA run
run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
        common_all_inputs=i_simple, unique_pt_inputs=i_arm,
        init_event_list=init_event_list, evt_react_list=evt_react_list,
        util_ongoing_list=util_ongoing, cost_ongoing_list=cost_ongoing, ipd=1,
        sensitivity_inputs=i_sensitivity, sensitivity_names=c("DSA_min","DSA_max"),
        sensitivity_bool=TRUE, n_sensitivity=length(l_inputs[[1]]),
        input_out=unlist(l_inputs[["parameter_name"]]))

# Scenario run
run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
        common_all_inputs=i_simple, unique_pt_inputs=i_arm,
        init_event_list=init_event_list, evt_react_list=evt_react_list,
        util_ongoing_list=util_ongoing, cost_ongoing_list=cost_ongoing, ipd=1,
        sensitivity_inputs=i_sensitivity, sensitivity_names=c("scenario_1","scenario_2"),
        sensitivity_bool=TRUE, n_sensitivity=1,
        input_out=unlist(l_inputs[["parameter_name"]]))
```

### Parameters spread across different levels

What happens when we have inputs not only at a single level (e.g.,
patient level), but also at simulation, patient, arm, etc levels? In
that case we can use
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md)
at each level, and the engine will automatically compute the right
offsets and accumulate the total `n_sensitivity` across all blocks — no
manual indicator tracking is needed.

For example, if we have 7 parameters at simulation level and 2 at
patient level, the engine knows the patient-level block starts at index
8 and applies the correct `n_sens_before` offset internally. The
[`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md)
helper illustrates how the indexing works:

``` r


#let's set the index to 8th (so it would correspond to the patient level indicator)
create_indicators(8,20,c(1,1,1,1,1,1,1),0) #all 0s
#> [1] 0 0 0 0 0 0 0

create_indicators(8,20,c(1,1,1),7) #first index is 1! because we know we are at index 8, and we have already gone through 7 iterations
#> [1] 1 0 0
```

Let’s simply add a few parameters that are also now at the patient
level. Note that as we are not using these parameters in the model, they
will have no impact but we can still see they are really being varied.
With
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md)
at each level, `n_sensitivity` and `sensitivity_names` are auto-detected
and the DSA offsets are computed automatically.

``` r


l_inputs_pat <- list(parameter_name = list("age","sex"),
                 base_value = list(60,1),
                 PSA_dist = list("rnorm","rbinom"),
                 a=list(60,1),
                 b=list(10,0.5),
                 n=as.list(rep(1,2)),
                 DSA_min = list(30,0),
                 DSA_max = list(80,1),
                 scenario_1=list(55,1),
                 scenario_2=list(45,0),
                 psa_indicators = as.list(rep(1,2))
                 )

i_simple <- input_block(
  base           = l_inputs[["base_value"]],
  psa            = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
  sens           = l_inputs,
  names_out      = l_inputs[["parameter_name"]],
  psa_indicators = l_inputs[["psa_indicators"]],
  dsa_names      = c("DSA_min", "DSA_max"))

i_pat <- input_block(
  base           = l_inputs_pat[["base_value"]],
  psa            = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]], l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
  sens           = l_inputs_pat,
  names_out      = l_inputs_pat[["parameter_name"]],
  psa_indicators = l_inputs_pat[["psa_indicators"]],
  dsa_names      = c("DSA_min", "DSA_max"))

results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  common_pt_inputs  = i_pat,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1,
  sensitivity_bool = TRUE,
  input_out = c(unlist(l_inputs[["parameter_name"]]), unlist(l_inputs_pat[["parameter_name"]]))
)
#> n_sensitivity auto-detected as 9 from input_block metadata
#> sensitivity_names auto-set to c("DSA_min", "DSA_max") from input_block metadata
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 1: 0.13s
#> Analysis number: 2
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 2: 0.13s
#> Analysis number: 3
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 3: 0.13s
#> Analysis number: 4
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 4: 0.13s
#> Analysis number: 5
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 5: 0.13s
#> Analysis number: 6
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 6: 0.13s
#> Analysis number: 7
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 7: 0.14s
#> Analysis number: 8
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 8: 0.13s
#> Analysis number: 9
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 9: 0.14s
#> Analysis number: 10
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 10: 0.13s
#> Analysis number: 11
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 11: 0.13s
#> Analysis number: 12
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 12: 0.14s
#> Analysis number: 13
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 13: 0.13s
#> Analysis number: 14
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 14: 0.13s
#> Analysis number: 15
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 15: 0.13s
#> Analysis number: 16
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 16: 0.13s
#> Analysis number: 17
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 17: 0.14s
#> Analysis number: 18
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 18: 0.13s
#> Total time to run: 2.38s
#> Simulation finalized;

summary_results_sens(results)
#>          arm analysis analysis_name     variable                value
#>       <char>    <int>        <char>       <fctr>               <char>
#>    1:    int        1       DSA_min        costs 7,768 (7,768; 7,768)
#>    2:  noint        1       DSA_min        costs 5,826 (5,826; 5,826)
#>    3:    int        2       DSA_min        costs 7,768 (7,768; 7,768)
#>    4:  noint        2       DSA_min        costs 5,826 (5,826; 5,826)
#>    5:    int        3       DSA_min        costs 3,884 (3,884; 3,884)
#>   ---                                                                
#> 1580:  noint       16       DSA_max dutil.sicker             0 (0; 0)
#> 1581:    int       17       DSA_max dutil.sicker             0 (0; 0)
#> 1582:  noint       17       DSA_max dutil.sicker             0 (0; 0)
#> 1583:    int       18       DSA_max dutil.sicker             0 (0; 0)
#> 1584:  noint       18       DSA_max dutil.sicker             0 (0; 0)


data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))


data_sensitivity |> group_by(sensitivity) |> summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int","age","sex"),mean)  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| sensitivity | util.sick | util.sicker | cost.sick | cost.sicker | cost.int | coef_noint | HR_int | age | sex |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.6 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 2 | 0.8 | 0.3 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 3 | 0.8 | 0.5 | 1000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 4 | 0.8 | 0.5 | 3000 | 5000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 5 | 0.8 | 0.5 | 3000 | 7000 | 800 | -1.609 | 0.8 | 60 | 1 |
| 6 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -2.303 | 0.8 | 60 | 1 |
| 7 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.5 | 60 | 1 |
| 8 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 30 | 1 |
| 9 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 0 |
| 10 | 0.9 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 11 | 0.8 | 0.7 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 12 | 0.8 | 0.5 | 5000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 13 | 0.8 | 0.5 | 3000 | 9000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 14 | 0.8 | 0.5 | 3000 | 7000 | 2000 | -1.609 | 0.8 | 60 | 1 |
| 15 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -0.916 | 0.8 | 60 | 1 |
| 16 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.9 | 60 | 1 |
| 17 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 80 | 1 |
| 18 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |

**Alternative: legacy approach using `pick_val_v()`**

``` r

# Manual tracking of n_sens_before and total n_sensitivity needed
i_sensitivity <- add_item(
  iterator_sensitivity = sens_iterator(sens, n_sensitivity)
) |>
  add_item(
    indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min", "DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                        rep(1, length(l_inputs[[1]])))
    } else { rep(1, length(l_inputs[[1]])) }
  ) |>
  add_item(
    indicators_pat = if (sensitivity_bool & sens_name_used %in% c("DSA_min", "DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                        rep(1, length(l_inputs_pat[[1]])), length(l_inputs[[1]]))
    } else { rep(1, length(l_inputs_pat[[1]])) }
  )

i_simple <- add_item() |>
  add_item(
    pick_val_v(base = l_inputs[["base_value"]],
               psa = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
               sens = l_inputs[[sens_name_used]], psa_ind = psa_bool, sens_ind = sensitivity_bool,
               indicator = indicators, names_out = l_inputs[["parameter_name"]],
               indicator_psa = l_inputs[["psa_indicators"]]))

i_pat <- add_item() |>
  add_item(
    pick_val_v(base = l_inputs_pat[["base_value"]],
               psa = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]], l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
               sens = l_inputs_pat[[sens_name_used]], psa_ind = psa_bool, sens_ind = sensitivity_bool,
               indicator = indicators_pat, names_out = l_inputs_pat[["parameter_name"]],
               indicator_psa = l_inputs_pat[["psa_indicators"]]))

run_sim(npats=5, n_sim=2, psa_bool=FALSE, arm_list=c("int","noint"),
        common_all_inputs=i_simple, unique_pt_inputs=i_arm, common_pt_inputs=i_pat,
        init_event_list=init_event_list, evt_react_list=evt_react_list,
        util_ongoing_list=util_ongoing, cost_ongoing_list=cost_ongoing, ipd=1,
        sensitivity_inputs=i_sensitivity, sensitivity_names=c("DSA_min","DSA_max"),
        sensitivity_bool=TRUE, n_sensitivity=length(l_inputs[[1]])+length(l_inputs_pat[[1]]),
        input_out=c(unlist(l_inputs[["parameter_name"]]),unlist(l_inputs_pat[["parameter_name"]])))
```

### Multiple parameters covaried

WARDEN also allows to have parameters being changed together in the DSA,
so if the programmer believes that e.g, age and sex should take their
min and max value together, that can be done by switching from a 0-1
vector to a vector which takes integer values to reflect the DSA
scenario number. See below this applied, and we also assume that the
utilities and costs are varied together. In this case, the indicators
are simplified, as we only need 1) the correct index and 2) the
`dsa_indicators` index to understand which parameters are covaried.
`n_sensitivity` is now auto-detected from the
[`input_block()`](https://jsanchezalv.github.io/WARDEN/reference/input_block.md)
metadata — no manual specification needed.

``` r


l_inputs <- list(parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
                 base_value = list(0.8,0.5,3000,7000,1000,log(0.2),0.8),
                 PSA_dist = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
                 a=list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)),
                 b=lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)), function(x) abs(x/5)),
                 n=as.list(rep(1,7)),
                 DSA_min = list(0.6,0.3,1000,5000,800,log(0.1),0.5),
                 DSA_max = list(0.9,0.7,5000,9000,2000,log(0.4),0.9),
                 scenario_1=list(0.6,0.3,1000,5000,800,log(0.1),0.5),
                 scenario_2=list(0.9,0.7,5000,9000,2000,log(0.4),0.9),
                 psa_indicators = as.list(c(rep(1,4),rep(0,3))),
                 dsa_indicators = as.list(c(1,1,2,2,2,3,4))  #covary utilities together, costs together
                 )

l_inputs_pat <- list(parameter_name = list("age","sex"),
                 base_value = list(60,1),
                 PSA_dist = list("rnorm","rbinom"),
                 a=list(60,1),
                 b=list(10,0.5),
                 n=as.list(rep(1,2)),
                 DSA_min = list(30,0),
                 DSA_max = list(80,1),
                 scenario_1=list(55,1),
                 scenario_2=list(45,0),
                 psa_indicators = as.list(rep(1,2)),
                 dsa_indicators = list(5,5) #will be 5th analysis done
                 )

i_simple <- input_block(
  base                  = l_inputs[["base_value"]],
  psa                   = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
  sens                  = l_inputs,
  names_out             = l_inputs[["parameter_name"]],
  psa_indicators        = l_inputs[["psa_indicators"]],
  dsa_names             = c("DSA_min", "DSA_max"),
  sens_indicators       = l_inputs[["dsa_indicators"]],
  indicator_sens_binary = FALSE,
  distributions         = l_inputs[["PSA_dist"]],
  covariances           = l_inputs[["b"]])

i_pat <- input_block(
  base                  = l_inputs_pat[["base_value"]],
  psa                   = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]], l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
  sens                  = l_inputs_pat,
  names_out             = l_inputs_pat[["parameter_name"]],
  psa_indicators        = l_inputs_pat[["psa_indicators"]],
  dsa_names             = c("DSA_min", "DSA_max"),
  sens_indicators       = l_inputs_pat[["dsa_indicators"]],
  indicator_sens_binary = FALSE,
  distributions         = l_inputs_pat[["PSA_dist"]],
  covariances           = l_inputs_pat[["b"]])


results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  common_pt_inputs = i_pat,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1,
  sensitivity_bool = TRUE,
  input_out = c(unlist(l_inputs[["parameter_name"]]), unlist(l_inputs_pat[["parameter_name"]]))
)
#> n_sensitivity auto-detected as 5 from input_block metadata
#> sensitivity_names auto-set to c("DSA_min", "DSA_max") from input_block metadata
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 1: 0.13s
#> Analysis number: 2
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 2: 0.13s
#> Analysis number: 3
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 3: 0.14s
#> Analysis number: 4
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 4: 0.13s
#> Analysis number: 5
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 5: 0.13s
#> Analysis number: 6
#> Simulation number: 1
#> Time to run simulation 1: 0.1s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 6: 0.16s
#> Analysis number: 7
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 7: 0.13s
#> Analysis number: 8
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 8: 0.13s
#> Analysis number: 9
#> Simulation number: 1
#> Time to run simulation 1: 0.08s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 9: 0.14s
#> Analysis number: 10
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 10: 0.13s
#> Total time to run: 1.37s
#> Simulation finalized;

summary_results_sens(results)
#>         arm analysis analysis_name     variable                value
#>      <char>    <int>        <char>       <fctr>               <char>
#>   1:    int        1       DSA_min        costs 7,768 (7,768; 7,768)
#>   2:  noint        1       DSA_min        costs 5,826 (5,826; 5,826)
#>   3:    int        2       DSA_min        costs 3,496 (3,496; 3,496)
#>   4:  noint        2       DSA_min        costs 1,942 (1,942; 1,942)
#>   5:    int        3       DSA_min        costs 7,768 (7,768; 7,768)
#>  ---                                                                
#> 876:  noint        8       DSA_max dutil.sicker             0 (0; 0)
#> 877:    int        9       DSA_max dutil.sicker             0 (0; 0)
#> 878:  noint        9       DSA_max dutil.sicker             0 (0; 0)
#> 879:    int       10       DSA_max dutil.sicker             0 (0; 0)
#> 880:  noint       10       DSA_max dutil.sicker             0 (0; 0)


data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))


data_sensitivity |> group_by(sensitivity) |> summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int","age","sex"),mean)  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| sensitivity | util.sick | util.sicker | cost.sick | cost.sicker | cost.int | coef_noint | HR_int | age | sex |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.6 | 0.3 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 2 | 0.8 | 0.5 | 1000 | 5000 | 800 | -1.609 | 0.8 | 60 | 1 |
| 3 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -2.303 | 0.8 | 60 | 1 |
| 4 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.5 | 60 | 1 |
| 5 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 30 | 0 |
| 6 | 0.9 | 0.7 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 |
| 7 | 0.8 | 0.5 | 5000 | 9000 | 2000 | -1.609 | 0.8 | 60 | 1 |
| 8 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -0.916 | 0.8 | 60 | 1 |
| 9 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.9 | 60 | 1 |
| 10 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 80 | 1 |

Old approach using `pick_val_v()` directly

``` r


i_sensitivity <- add_item(
  iterator_sensitivity = sens_iterator(sens, n_sensitivity)
)

i_simple <- add_item() |>
  add_item(
    pick_val_v(
      base = l_inputs[["base_value"]],
      psa = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
      sens          = l_inputs[[sens_name_used]],
      psa_ind       = psa_bool,
      sens_ind      = sensitivity_bool,
      indicator     = l_inputs[["dsa_indicators"]],
      sens_iterator = iterator_sensitivity,
      indicator_sens_binary = FALSE,
      names_out     = l_inputs[["parameter_name"]],
      indicator_psa = l_inputs[["psa_indicators"]],
      distributions = l_inputs[["PSA_dist"]],
      covariances   = l_inputs[["b"]]
    )
  )

i_pat <- add_item() |>
  add_item(
    pick_val_v(
      base = l_inputs_pat[["base_value"]],
      psa = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]], l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
      sens          = l_inputs_pat[[sens_name_used]],
      psa_ind       = psa_bool,
      sens_ind      = sensitivity_bool,
      indicator     = l_inputs_pat[["dsa_indicators"]],
      sens_iterator = iterator_sensitivity,
      indicator_sens_binary = FALSE,
      names_out     = l_inputs_pat[["parameter_name"]],
      indicator_psa = l_inputs_pat[["psa_indicators"]],
      distributions = l_inputs_pat[["PSA_dist"]],
      covariances   = l_inputs_pat[["b"]]
    )
  )

results <- run_sim(
  npats=5, n_sim=2, psa_bool = FALSE, arm_list = c("int", "noint"),
  common_all_inputs = i_simple, unique_pt_inputs = i_arm, common_pt_inputs = i_pat,
  init_event_list = init_event_list, evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing, cost_ongoing_list = cost_ongoing, ipd = 1,
  sensitivity_inputs = i_sensitivity,
  sensitivity_names = c("DSA_min", "DSA_max"),
  sensitivity_bool = TRUE,
  n_sensitivity = length(unique(unlist(l_inputs[["dsa_indicators"]]))) +
                  length(unique(unlist(l_inputs_pat[["dsa_indicators"]]))), #5 parameters
  input_out = c(unlist(l_inputs[["parameter_name"]]), unlist(l_inputs_pat[["parameter_name"]]))
)
```

### Parameters that are vectors

Now we are ready to handle all cases. The only approach not shown yet is
when some parameters are vectors instead of length 1. In this case, the
approach is the same as the previous one, as those parameters would be
varied together. We add to `l_inputs_pat` a new parameter of length 2
that will be a multivariate normal. Because `l_inputs_pat` changes, we
rebuild `i_pat` accordingly — `i_simple` and
[`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
are otherwise unchanged. Below it can be seen for a DSA, but by now it
should be straightforward to switch to probabilistic DSA, deterministic
case, or standard PSA.

``` r


l_inputs_pat <- list(parameter_name = list("age","sex", "v_state"),
                 base_value = list(60,1, c(10,20)),
                 PSA_dist = list("rnorm","rbinom","mvrnorm"),
                 a=list(60,1,c(10,20)),
                 b=list(10,0.5,matrix(c(2,1,4,1),2,2)),
                 n=as.list(rep(1,3)),
                 DSA_min = list(30,0, c(5,10)),
                 DSA_max = list(80,1, c(15,25)),
                 scenario_1=list(55,1, c(12,21)),
                 scenario_2=list(45,0, c(16,10)),
                 psa_indicators = list(1,1,c(1,0)), #we vary the first one but not the second one in PSA!
                 dsa_indicators = list(5,5,c(6,6))
                 )

i_pat <- input_block(
  base                  = l_inputs_pat[["base_value"]],
  psa                   = pick_psa(l_inputs_pat[["PSA_dist"]], l_inputs_pat[["n"]], l_inputs_pat[["a"]], l_inputs_pat[["b"]]),
  sens                  = l_inputs_pat,
  names_out             = l_inputs_pat[["parameter_name"]],
  psa_indicators        = l_inputs_pat[["psa_indicators"]],
  dsa_names             = c("DSA_min", "DSA_max"),
  sens_indicators       = l_inputs_pat[["dsa_indicators"]],
  indicator_sens_binary = FALSE,
  distributions         = l_inputs_pat[["PSA_dist"]],
  covariances           = l_inputs_pat[["b"]])

results <- run_sim(
  npats=5,
  n_sim=2,
  psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = i_simple,
  unique_pt_inputs  = i_arm,
  common_pt_inputs = i_pat,
  init_event_list = init_event_list,
  evt_react_list = evt_react_list,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1,
  sensitivity_bool = TRUE,
  input_out = c(unlist(l_inputs[["parameter_name"]]), unlist(l_inputs_pat[["parameter_name"]]))
)
#> n_sensitivity auto-detected as 6 from input_block metadata
#> sensitivity_names auto-set to c("DSA_min", "DSA_max") from input_block metadata
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 1: 0.13s
#> Analysis number: 2
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 2: 0.13s
#> Analysis number: 3
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 3: 0.14s
#> Analysis number: 4
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 4: 0.13s
#> Analysis number: 5
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 5: 0.14s
#> Analysis number: 6
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 6: 0.13s
#> Analysis number: 7
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 7: 0.13s
#> Analysis number: 8
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.07s
#> Time to run analysis 8: 0.14s
#> Analysis number: 9
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 9: 0.13s
#> Analysis number: 10
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 10: 0.13s
#> Analysis number: 11
#> Simulation number: 1
#> Time to run simulation 1: 0.07s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 11: 0.14s
#> Analysis number: 12
#> Simulation number: 1
#> Time to run simulation 1: 0.06s
#> Simulation number: 2
#> Time to run simulation 2: 0.06s
#> Time to run analysis 12: 0.13s
#> Total time to run: 1.61s
#> Simulation finalized;

summary_results_sens(results)
#>          arm analysis analysis_name     variable                value
#>       <char>    <int>        <char>       <fctr>               <char>
#>    1:    int        1       DSA_min        costs 7,768 (7,768; 7,768)
#>    2:  noint        1       DSA_min        costs 5,826 (5,826; 5,826)
#>    3:    int        2       DSA_min        costs 3,496 (3,496; 3,496)
#>    4:  noint        2       DSA_min        costs 1,942 (1,942; 1,942)
#>    5:    int        3       DSA_min        costs 7,768 (7,768; 7,768)
#>   ---                                                                
#> 1052:  noint       10       DSA_max dutil.sicker             0 (0; 0)
#> 1053:    int       11       DSA_max dutil.sicker             0 (0; 0)
#> 1054:  noint       11       DSA_max dutil.sicker             0 (0; 0)
#> 1055:    int       12       DSA_max dutil.sicker             0 (0; 0)
#> 1056:  noint       12       DSA_max dutil.sicker             0 (0; 0)


data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))

#v_state is length > 1, so it does not appear in merged_df as it would break the data. We can get it manually and print it
v_state_avg <- tibble( v_state =
    unlist(
      lapply(
        map_depth(results,2, "v_state"), function(x) {
            vecs <- list(x[[1]]$int, x[[1]]$noint, x[[2]]$int, x[[2]]$noint)
            paste(colMeans(do.call(rbind, vecs)), collapse=", ")
        }
      )
    )
  )

data_sensitivity |> group_by(sensitivity) |> summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int","age","sex"),mean) |>
  bind_cols(v_state_avg)  |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

| sensitivity | util.sick | util.sicker | cost.sick | cost.sicker | cost.int | coef_noint | HR_int | age | sex | v_state |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| 1 | 0.6 | 0.3 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 | 10, 20 |
| 2 | 0.8 | 0.5 | 1000 | 5000 | 800 | -1.609 | 0.8 | 60 | 1 | 10, 20 |
| 3 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -2.303 | 0.8 | 60 | 1 | 10, 20 |
| 4 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.5 | 60 | 1 | 10, 20 |
| 5 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 30 | 0 | 10, 20 |
| 6 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 | 5, 10 |
| 7 | 0.9 | 0.7 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 | 10, 20 |
| 8 | 0.8 | 0.5 | 5000 | 9000 | 2000 | -1.609 | 0.8 | 60 | 1 | 10, 20 |
| 9 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -0.916 | 0.8 | 60 | 1 | 10, 20 |
| 10 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.9 | 60 | 1 | 10, 20 |
| 11 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 80 | 1 | 10, 20 |
| 12 | 0.8 | 0.5 | 3000 | 7000 | 1000 | -1.609 | 0.8 | 60 | 1 | 15, 25 |
