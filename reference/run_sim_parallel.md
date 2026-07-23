# Run simulations in parallel mode (at the simulation level)

Run simulations in parallel mode (at the simulation level)

## Usage

``` r
run_sim_parallel(
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
  ncores = 1,
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

- ncores:

  The number of cores to use for parallel computing

- input_out:

  A vector of variables to be returned in the output data frame

- ipd:

  Integer taking value 0 if no IPD data returned, 1 for full IPD data
  returned, and 2 IPD data but aggregating events

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
  it is updated. Recommended for resource-constrained models where queue
  waits make future event times unpredictable. With TRUE, ongoing values
  are computed at each event for the interval just completed, avoiding
  prospective errors from
  [`adj_val()`](https://jsanchezalv.github.io/WARDEN/reference/adj_val.md).

- continue_on_error:

  If TRUE, on error at patient stage will attempt to continue to the
  next simulation (only works if n_sim and/or n_sensitivity are \> 1,
  not at the patient level)

- seed:

  Starting seed to be used for the whole analysis. If null, it's set to
  1 by default.

## Value

A list of lists with the analysis results

## Details

This function is slightly different from `run_sim`. `run_sim` allows to
run single-core. `run_sim_parallel` allows to use multiple-core at the
simulation level, making it more efficient for a large number of
simulations relative to `run_sim` (e.g., for PSA).

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
that if ncores \> 1, then results per simulation will only be exactly
replicable if using run_sim_parallel (as seeds are automatically
transformed to be seven integer seeds -i.e, L'Ecuyer-CMRG seeds-) Note
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

If `continue_on_error` is set to FALSE, it will only export analysis
level inputs due to the parallel engine (use single-engine for those
inputs) `continue_on_error` will skip the current simulation (so it
won't continue for the rest of patient-arms) if TRUE. Note that this
will make the progress bar not correct, as a set of patients that were
expected to be run is not.

## Examples

``` r
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
             input = {}) |>
  add_reactevt(name_evt = "sicker",
               input = {
                 q_default <- util.sicker
                 c_default <- cost.sicker + if(arm=="int"){cost.int}else{0}
                 fl.sick <- 0
               }) |>
  add_reactevt(name_evt = "death",
               input = {
                 q_default <- 0
                 c_default <- 0
                 curtime <- Inf
               }) 
               
util_ongoing <- "q_default"
cost_ongoing <- "c_default"
                          

run_sim_parallel(arm_list=c("int","noint"),
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
ipd = 1,
ncores = 1)
#> Analysis number: 1
#> Simulation number: 1
#> Time to run simulation 1: 0.05s
#> Time to run analysis 1: 0.76s
#> Total time to run: 0.76s
#> Simulation finalized; 
#> WARDEN Simulation Results
#> -------------------------
#>   Analyses (sensitivities): 1 
#>   Simulations per analysis: 1 
#>   Arms: int, noint 
#> 
#> Access patterns:
#>   Single sim:  results[[sens]][[sim]]   (e.g., results[[1]][[1]])
#>   Sim list:    results[[sens]]           (e.g., results[[1]])
#>   Full object: results
#> 
#> Summary functions:
#>   summary_results_det(results, sens=1, sim=1)   # deterministic
#>   summary_results_sim(results, sens=1)          # PSA
#>   summary_results_sens(results)                 # sensitivity
```
