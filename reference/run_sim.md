# Run the simulation

Run the simulation

## Usage

``` r
run_sim(
  arm_list = c("int", "noint"),
  sensitivity_inputs = NULL,
  common_all_inputs = NULL,
  common_pt_inputs = NULL,
  unique_pt_inputs = NULL,
  init_event_list = NULL,
  evt_react_list = evt_react_list,
  util_ongoing_list = NULL,
  util_instant_list = NULL,
  util_cycle_list = NULL,
  cost_ongoing_list = NULL,
  cost_instant_list = NULL,
  cost_cycle_list = NULL,
  other_ongoing_list = NULL,
  other_instant_list = NULL,
  npats = 500,
  n_sim = 1,
  psa_bool = NULL,
  sensitivity_bool = FALSE,
  sensitivity_names = NULL,
  n_sensitivity = 1,
  input_out = character(),
  ipd = 1,
  constrained = FALSE,
  timed_freq = NULL,
  debug = FALSE,
  accum_backwards = FALSE,
  continue_on_error = FALSE,
  seed = NULL
)
```

## Arguments

- arm_list:

  A vector of the names of the interventions evaluated in the simulation

- sensitivity_inputs:

  A list of sensitivity inputs that do not change within a sensitivity
  in a similar fashion to common_all_inputs, etc

- common_all_inputs:

  A list of inputs common across patients that do not change within a
  simulation

- common_pt_inputs:

  A list of inputs that change across patients but are not affected by
  the intervention

- unique_pt_inputs:

  A list of inputs that change across each intervention

- init_event_list:

  A list of initial events and event times. If no initial events are
  given, a "Start" event at time 0 is created automatically

- evt_react_list:

  A list of event reactions

- util_ongoing_list:

  Vector of QALY named variables that are accrued at an ongoing basis
  (discounted using drq)

- util_instant_list:

  Vector of QALY named variables that are accrued instantaneously at an
  event (discounted using drq)

- util_cycle_list:

  Vector of QALY named variables that are accrued in cycles (discounted
  using drq)

- cost_ongoing_list:

  Vector of cost named variables that are accrued at an ongoing basis
  (discounted using drc)

- cost_instant_list:

  Vector of cost named variables that are accrued instantaneously at an
  event (discounted using drc)

- cost_cycle_list:

  Vector of cost named variables that are accrued in cycles (discounted
  using drc)

- other_ongoing_list:

  Vector of other named variables that are accrued at an ongoing basis
  (discounted using drq)

- other_instant_list:

  Vector of other named variables that are accrued instantaneously at an
  event (discounted using drq)

- npats:

  The number of patients to be simulated (it will simulate npats \*
  length(arm_list))

- n_sim:

  The number of simulations to run per sensitivity

- psa_bool:

  A boolean to determine if PSA should be conducted. If n_sim \> 1 and
  psa_bool = FALSE, the differences between simulations will be due to
  sampling

- sensitivity_bool:

  A boolean to determine if Scenarios/DSA should be conducted.

- sensitivity_names:

  A vector of scenario/DSA names that can be used to select the right
  sensitivity (e.g., c("Scenario_1", "Scenario_2")). The parameter
  "sens_name_used" is created from it which corresponds to the one being
  used for each iteration.

- n_sensitivity:

  Number of sensitivity analysis (DSA or Scenarios) to run. It will be
  interacted with sensitivity_names argument if not null
  (n_sensitivityitivity = n_sensitivity \* length(sensitivity_names)).
  For DSA, it should be as many parameters as there are. For scenario,
  it should be 1.

- input_out:

  A vector of variables to be returned in the output data frame

- ipd:

  Integer taking value 1 for full IPD data returned, and 2 IPD data but
  aggregating events (returning last value for numeric/character/factor
  variables. For other objects (e.g., matrices), the IPD will still be
  returned as the aggregation rule is not clear). Other values mean no
  IPD data returned (removes non-numerical or length\>1 items)

- constrained:

  Boolean, FALSE by default, which runs the simulation with patients not
  interacting with each other, TRUE if resources are shared within an
  arm (allows constrained resources)

- timed_freq:

  If NULL, it does not produce any timed outputs. Otherwise should be a
  number (e.g., every 1 year)

- debug:

  If TRUE, will generate a log file

- accum_backwards:

  If TRUE, the ongoing accumulators will count backwards (i.e., the
  current value is applied until the previous update). If FALSE, the
  current value is applied between the current event and the next time
  it is updated.

- continue_on_error:

  If TRUE, on error it will attempt to continue by skipping the current
  simulation

- seed:

  Starting seed to be used for the whole analysis. If null, it's set to
  1 by default.

## Value

A list of data frames with the simulation results

## Details

This function is slightly different from `run_sim_parallel`.
`run_sim_parallel` only runs multiple-core at the simulation level.
`run_sim` uses only-single core. `run_sim` can be more efficient if
using only one simulation (e.g., deterministic), while
`run_sim_parallel` will be more efficient if the number of simulations
is \>1 (e.g., PSA).

Event ties are processed in the order declared within the
`init_event_list` argument (`evts` argument within the first sublist of
that object). To do so, the program automatically adds a sequence from
to 0 to the (number of events - 1) times 1e-10 to add to the event times
when selecting the event with minimum time. This time has been selected
as it's relatively small yet not so small as to be ignored by which.min
(see .Machine for more details)

A list of protected objects that should not be used by the user as input
names or in the global environment to avoid the risk of overwriting them
is as follows: c("arm", "arm_list", "categories_for_export",
"cur_evtlist", "curtime", "evt", "i", "prevtime", "sens", "simulation",
"sens_name_used","list_env","uc_lists","npats","ipd").

The engine uses the L'Ecuyer-CMRG for the random number generator. Note
that the random seeds are set to be unique in their category (i.e., at
patient level, patient-arm level, etc.)

If no `drc` or `drq` parameters are passed within `sensitivity` or
`common_all` input lists, these are assigned a default value 0.03 for
discounting costs, QALYs and others.

Ongoing items will look backward to the last time updated when
performing the discounting and accumulation. This means that the user
does not necessarily need to keep updating the value, but only add it
when the value changes looking forward (e.g., o_q = utility at event 1,
at event 2 utility does not change, but at event 3 it does, so we want
to make sure to add o_q = utility at event 3 before updating utility.
The program will automatically look back until event 1). Note that in
previous versions of the package backward was the default, and now this
has switched to forward.

The requirement to use `modify_item` if using `accum_backwards = TRUE`,
is no longer the case thanks to a new method using active bindings, so
it can be used normally.

It is important to note that the QALYs and Costs (ongoing or instant or
per cycle) used should be of length 1. If they were of length \> 1, the
model would expand the data, so instead of having each event as a row,
the event would have N rows (equal to the length of the costs/qalys to
discount passed). This means more processing of the results data would
be needed in order for it to provide the correct results.

If the `cycle` lists are used, then it is expected the user will declare
as well the name of the variable pasted with `cycle_l` and
`cycle_starttime` (e.g., c_default_cycle_l and
c_default_cycle_starttime) to ensure the discounting can be computed
using cycles, with cycle_l being the cycle length, and cycle_starttime
being the starting time in which the variable started counting.
Optionally, `max_cycles` must also be added (if no maximum number of
cycles, it should be set equal to NA).

`debug = TRUE` will export a log file with the timestamp up the error in
the main working directory. Note that using this mode without
modify_item or modify_item_seq may lead to inaccuracies if assignments
are done in non-standard ways, as the AST may not catch all the relevant
assignments (e.g., an assigment like assign(paste("x\_",i),5) in a loop
will not be identified).

`continue_on_error` will skip the current simulation (so it won't
continue for the rest of patient-arms) if TRUE. Note that this will make
the progress bar not correct, as a set of patients that were expected to
be run is not.

## Examples

``` r
library(magrittr)
common_all_inputs <-add_item(
util.sick = 0.8,
util.sicker = 0.5,
cost.sick = 3000,
cost.sicker = 7000,
cost.int = 1000,
coef_noint = log(0.2),
HR_int = 0.8,
drc = 0.035, #different values than what's assumed by default
drq = 0.035,
random_seed_sicker_i = sample.int(100000,5,replace = FALSE)
)

common_pt_inputs <- add_item(death= max(0.0000001,rnorm(n=1, mean=12, sd=3))) 

unique_pt_inputs <- add_item(fl.sick = 1,
                             q_default = util.sick,
                             c_default = cost.sick + if(arm=="int"){cost.int}else{0}) 
                             
init_event_list <- 
add_tte(arm=c("noint","int"), evts = c("sick","sicker","death") ,input={
  sick <- 0
  sicker <- draw_tte(1,dist="exp",
   coef1=coef_noint, beta_tx = ifelse(arm=="int",HR_int,1),
    seed = random_seed_sicker_i[i])
  
})   

evt_react_list <-
add_reactevt(name_evt = "sick",
             input = {}) %>%
  add_reactevt(name_evt = "sicker",
               input = {
                 q_default <- util.sicker
                 c_default <- cost.sicker + if(arm=="int"){cost.int}else{0}
                 fl.sick <- 0 
               }) %>%
  add_reactevt(name_evt = "death",
               input = {
                 q_default <- 0
                 c_default <- 0
                 curtime <- Inf
               }) 
               
util_ongoing <- "q_default"
cost_ongoing <- "c_default"
                          

run_sim(arm_list=c("int","noint"),
common_all_inputs = common_all_inputs,
common_pt_inputs = common_pt_inputs,
unique_pt_inputs = unique_pt_inputs,
init_event_list = init_event_list,
evt_react_list = evt_react_list,
util_ongoing_list = util_ongoing,
cost_ongoing_list = cost_ongoing,
npats = 2,
n_sim = 1,
psa_bool = FALSE,
ipd = 1)
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.05s
#> Time to run analysis 1: 0.05s
#> Total time to run: 0.06s
#> Simulation finalized; 
#> [[1]]
#> [[1]][[1]]
#> [[1]][[1]]$sensitivity_name
#> [1] ""
#> 
#> [[1]][[1]]$arm_list
#> [1] "int"   "noint"
#> 
#> [[1]][[1]]$total_lys
#>      int    noint 
#> 9.046874 9.046874 
#> 
#> [[1]][[1]]$total_qalys
#>      int    noint 
#> 6.207438 6.181151 
#> 
#> [[1]][[1]]$total_costs
#>      int    noint 
#> 49921.64 41225.25 
#> 
#> [[1]][[1]]$total_lys_undisc
#>      int    noint 
#> 10.89866 10.89866 
#> 
#> [[1]][[1]]$total_qalys_undisc
#>      int    noint 
#> 7.501176 7.474146 
#> 
#> [[1]][[1]]$total_costs_undisc
#>      int    noint 
#> 59831.36 49293.10 
#> 
#> [[1]][[1]]$c_default
#>      int    noint 
#> 49921.64 41225.25 
#> 
#> [[1]][[1]]$c_default_undisc
#>      int    noint 
#> 59831.36 49293.10 
#> 
#> [[1]][[1]]$q_default
#>      int    noint 
#> 6.207438 6.181151 
#> 
#> [[1]][[1]]$q_default_undisc
#>      int    noint 
#> 7.501176 7.474146 
#> 
#> [[1]][[1]]$merged_df
#>     evtname    evttime  prevtime pat_id    arm total_lys total_qalys
#>      <char>      <num>     <num>  <int> <char>     <num>       <num>
#>  1:    sick  0.0000000 0.0000000      1    int 10.339480    8.271584
#>  2:   death 12.7779512 0.0000000      1    int 10.339480    8.271584
#>  3:    sick  0.0000000 0.0000000      2    int  7.754267    4.143293
#>  4:  sicker  0.9010175 0.0000000      2    int  7.754267    4.143293
#>  5:   death  9.0193725 0.9010175      2    int  7.754267    4.143293
#>  6:    sick  0.0000000 0.0000000      1  noint 10.339480    8.271584
#>  7:   death 12.7779512 0.0000000      1  noint 10.339480    8.271584
#>  8:    sick  0.0000000 0.0000000      2  noint  7.754267    4.090719
#>  9:  sicker  0.7208140 0.0000000      2  noint  7.754267    4.090719
#> 10:   death  9.0193725 0.7208140      2  noint  7.754267    4.090719
#>     total_costs total_costs_undisc total_qalys_undisc total_lys_undisc
#>           <num>              <num>              <num>            <num>
#>  1:    41357.92           51111.80          10.222361        12.777951
#>  2:    41357.92           51111.80          10.222361        12.777951
#>  3:    58485.35           68550.91           4.779991         9.019372
#>  4:    58485.35           68550.91           4.779991         9.019372
#>  5:    58485.35           68550.91           4.779991         9.019372
#>  6:    31018.44           38333.85          10.222361        12.777951
#>  7:    31018.44           38333.85          10.222361        12.777951
#>  8:    51432.07           60252.35           4.725930         9.019372
#>  9:    51432.07           60252.35           4.725930         9.019372
#> 10:    51432.07           60252.35           4.725930         9.019372
#>            lys     qalys     costs lys_undisc qalys_undisc costs_undisc
#>          <num>     <num>     <num>      <num>        <num>        <num>
#>  1: 10.3394801 8.2715841 41357.920 12.7779512   10.2223609    51111.805
#>  2:  0.0000000 0.0000000     0.000  0.0000000    0.0000000        0.000
#>  3:  0.8871965 0.7097572  3548.786  0.9010175    0.7208140     3604.070
#>  4:  6.8670706 3.4335353 54936.565  8.1183550    4.0591775    64946.840
#>  5:  0.0000000 0.0000000     0.000  0.0000000    0.0000000        0.000
#>  6: 10.3394801 8.2715841 31018.440 12.7779512   10.2223609    38333.854
#>  7:  0.0000000 0.0000000     0.000  0.0000000    0.0000000        0.000
#>  8:  0.7119504 0.5695603  2135.851  0.7208140    0.5766512     2162.442
#>  9:  7.0423168 3.5211584 49296.218  8.2985585    4.1492793    58089.910
#> 10:  0.0000000 0.0000000     0.000  0.0000000    0.0000000        0.000
#>     c_default q_default c_default_undisc q_default_undisc   nexttime simulation
#>         <num>     <num>            <num>            <num>      <num>      <int>
#>  1: 41357.920 8.2715841        51111.805       10.2223609 12.7779512          1
#>  2:     0.000 0.0000000            0.000        0.0000000 12.7779512          1
#>  3:  3548.786 0.7097572         3604.070        0.7208140  0.9010175          1
#>  4: 54936.565 3.4335353        64946.840        4.0591775  9.0193725          1
#>  5:     0.000 0.0000000            0.000        0.0000000  9.0193725          1
#>  6: 31018.440 8.2715841        38333.854       10.2223609 12.7779512          1
#>  7:     0.000 0.0000000            0.000        0.0000000 12.7779512          1
#>  8:  2135.851 0.5695603         2162.442        0.5766512  0.7208140          1
#>  9: 49296.218 3.5211584        58089.910        4.1492793  9.0193725          1
#> 10:     0.000 0.0000000            0.000        0.0000000  9.0193725          1
#>     sensitivity
#>           <int>
#>  1:           1
#>  2:           1
#>  3:           1
#>  4:           1
#>  5:           1
#>  6:           1
#>  7:           1
#>  8:           1
#>  9:           1
#> 10:           1
#> 
#> 
#> 
```
