# Validate model inputs before running

Performs static checks on the model inputs without running the
simulation. Catches common configuration errors early.

## Usage

``` r
validate_model(
  arm_list = c("int", "noint"),
  init_event_list = NULL,
  evt_react_list = NULL,
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
  seed = NULL
)
```

## Arguments

- arm_list:

  Character vector of intervention names.

- init_event_list:

  List of initial events and event times (output of
  [`add_tte()`](https://jsanchezalv.github.io/WARDEN/reference/add_tte.md)).

- evt_react_list:

  List of event reactions (output of
  [`add_reactevt()`](https://jsanchezalv.github.io/WARDEN/reference/add_reactevt.md)).

- util_ongoing_list:

  Character vector of ongoing QALY variable names.

- util_instant_list:

  Character vector of instant QALY variable names.

- util_cycle_list:

  Character vector of cycle QALY variable names.

- cost_ongoing_list:

  Character vector of ongoing cost variable names.

- cost_instant_list:

  Character vector of instant cost variable names.

- cost_cycle_list:

  Character vector of cycle cost variable names.

- other_ongoing_list:

  Character vector of other ongoing variable names.

- other_instant_list:

  Character vector of other instant variable names.

- npats:

  Integer. Number of patients per simulation.

- n_sim:

  Integer. Number of simulations.

- psa_bool:

  Logical. Whether PSA mode is intended.

- seed:

  Numeric. Random seed.

## Value

Invisibly `TRUE` if all checks pass. Prints a summary of checks.
