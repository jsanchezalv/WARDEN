# Package index

## Main Modelling Functions

### Core Functions

- [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
  : Run the simulation
- [`run_sim_parallel()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim_parallel.md)
  : Run simulations in parallel mode (at the simulation level)
- [`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md)
  : Define or append model inputs
- [`add_reactevt()`](https://jsanchezalv.github.io/WARDEN/reference/add_reactevt.md)
  : Define the modifications to other events, costs, utilities, or other
  items affected by the occurrence of the event
- [`add_tte()`](https://jsanchezalv.github.io/WARDEN/reference/add_tte.md)
  : Define events and the initial event time
- [`new_event()`](https://jsanchezalv.github.io/WARDEN/reference/new_event.md)
  : Add Events to the Queue for a Patient
- [`modify_event()`](https://jsanchezalv.github.io/WARDEN/reference/modify_event.md)
  : Modify or Add Events for a Patient
- [`remove_event()`](https://jsanchezalv.github.io/WARDEN/reference/remove_event.md)
  : Remove Events for a Patient
- [`get_event()`](https://jsanchezalv.github.io/WARDEN/reference/get_event.md)
  : Get a specific event time
- [`has_event()`](https://jsanchezalv.github.io/WARDEN/reference/has_event.md)
  : Check if a Patient Has a Specific Event
- [`next_event()`](https://jsanchezalv.github.io/WARDEN/reference/next_event.md)
  : Get the Next Event(s) in the Queue
- [`next_event_pt()`](https://jsanchezalv.github.io/WARDEN/reference/next_event_pt.md)
  : Get the Next Event(s) in the Queue for a specific patient
- [`queue_empty()`](https://jsanchezalv.github.io/WARDEN/reference/queue_empty.md)
  : Check if the Event Queue is Empty
- [`queue_size()`](https://jsanchezalv.github.io/WARDEN/reference/queue_size.md)
  : Get the Size of the Event Queue
- [`queue_create()`](https://jsanchezalv.github.io/WARDEN/reference/queue_create.md)
  : Create a New Event Queue
- [`pop_event()`](https://jsanchezalv.github.io/WARDEN/reference/pop_event.md)
  : Remove the Next Event from the Queue
- [`pop_and_return_event()`](https://jsanchezalv.github.io/WARDEN/reference/pop_and_return_event.md)
  : Pop and Return the Next Event

### Resource Constrained Specific Functions

- [`resource_discrete()`](https://jsanchezalv.github.io/WARDEN/reference/resource_discrete.md)
  : Create a Discrete Resource
- [`shared_input()`](https://jsanchezalv.github.io/WARDEN/reference/shared_input.md)
  : Shared input object
- [`discrete_resource_clone()`](https://jsanchezalv.github.io/WARDEN/reference/discrete_resource_clone.md)
  : Clone independent discrete resources

### Auxiliary Functions

- [`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md)
  : Creates a vector of indicators (0 and 1) for sensitivity/DSA
  analysis

- [`sens_iterator()`](https://jsanchezalv.github.io/WARDEN/reference/sens_iterator.md)
  : Create an iterator based on sens of the current iteration within a
  scenario (DSA)

- [`pick_psa()`](https://jsanchezalv.github.io/WARDEN/reference/pick_psa.md)
  :

  Helper function to create a list with random draws or whenever a
  series of functions needs to be called. Can be implemented within
  `pick_val_v`.

- [`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
  : Select which values should be applied in the corresponding loop for
  several values (vector or list).

- [`replicate_profiles()`](https://jsanchezalv.github.io/WARDEN/reference/replicate_profiles.md)
  : Replicate profiles data.frame

### Summary Functions

- [`ceac_des()`](https://jsanchezalv.github.io/WARDEN/reference/ceac_des.md)
  : Calculate the cost-effectiveness acceptability curve (CEAC) for a
  DES model with a PSA result
- [`evpi_des()`](https://jsanchezalv.github.io/WARDEN/reference/evpi_des.md)
  : Calculate the Expected Value of Perfect Information (EVPI) for a DES
  model with a PSA result
- [`summary_results_det()`](https://jsanchezalv.github.io/WARDEN/reference/summary_results_det.md)
  : Deterministic results for a specific treatment
- [`summary_results_sens()`](https://jsanchezalv.github.io/WARDEN/reference/summary_results_sens.md)
  : Summary of sensitivity outputs for a treatment
- [`summary_results_sim()`](https://jsanchezalv.github.io/WARDEN/reference/summary_results_sim.md)
  : Summary of PSA outputs for a treatment

## Statistics Functions

### Auxiliary Functions

- [`draw_tte()`](https://jsanchezalv.github.io/WARDEN/reference/draw_tte.md)
  : Draw a time to event from a list of parametric survival functions
- [`random_stream()`](https://jsanchezalv.github.io/WARDEN/reference/random_stream.md)
  : Creates an environment (similar to R6 class) of random uniform
  numbers to be drawn from
- [`luck_adj()`](https://jsanchezalv.github.io/WARDEN/reference/luck_adj.md)
  : Perform luck adjustment
- [`qtimecov()`](https://jsanchezalv.github.io/WARDEN/reference/qtimecov.md)
  : Draw Time-to-Event with Time-Dependent Covariates and Luck
  Adjustment
- [`adj_val()`](https://jsanchezalv.github.io/WARDEN/reference/adj_val.md)
  : Adjusted Value Calculation

### Distributions Functions

- [`cond_dirichlet()`](https://jsanchezalv.github.io/WARDEN/reference/cond_dirichlet.md)
  : Calculate conditional dirichlet values
- [`cond_mvn()`](https://jsanchezalv.github.io/WARDEN/reference/cond_mvn.md)
  : Calculate conditional multivariate normal values
- [`rgamma_mse()`](https://jsanchezalv.github.io/WARDEN/reference/rgamma_mse.md)
  : Draw from a gamma distribution based on mean and se
- [`qgamma_mse()`](https://jsanchezalv.github.io/WARDEN/reference/qgamma_mse.md)
  : Use quantiles from a gamma distribution based on mean and se
- [`rpoisgamma()`](https://jsanchezalv.github.io/WARDEN/reference/rpoisgamma.md)
  : Draw time to event (tte) from a Poisson or Poisson-Gamma (PG)
  Mixture/Negative Binomial (NB) Process
- [`rpoisgamma_rcpp()`](https://jsanchezalv.github.io/WARDEN/reference/rpoisgamma_rcpp.md)
  : Draw time to event (tte) from a Poisson or Poisson-Gamma (PG)
  Mixture/Negative Binomial (NB) Process using C++
- [`qbeta_mse()`](https://jsanchezalv.github.io/WARDEN/reference/qbeta_mse.md)
  : Draw from a beta distribution based on mean and se (quantile)
- [`qcond_exp()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_exp.md)
  : Conditional quantile function for exponential distribution
- [`qcond_gamma()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_gamma.md)
  : Conditional quantile function for gamma distribution
- [`qcond_gompertz()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_gompertz.md)
  : Quantile function for conditional Gompertz distribution (lower bound
  only)
- [`qcond_llogis()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_llogis.md)
  : Conditional quantile function for loglogistic distribution
- [`qcond_lnorm()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_lnorm.md)
  : Conditional quantile function for lognormal distribution
- [`qcond_norm()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_norm.md)
  : Conditional quantile function for normal distribution
- [`qcond_weibull()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_weibull.md)
  : Conditional quantile function for weibull distribution
- [`qcond_weibullPH()`](https://jsanchezalv.github.io/WARDEN/reference/qcond_weibullPH.md)
  : Conditional quantile function for WeibullPH (flexsurv)
- [`rbeta_mse()`](https://jsanchezalv.github.io/WARDEN/reference/rbeta_mse.md)
  : Draw from a beta distribution based on mean and se
- [`rcond_gompertz()`](https://jsanchezalv.github.io/WARDEN/reference/rcond_gompertz.md)
  : Draw from a conditional Gompertz distribution (lower bound only)
- [`rcond_gompertz_lu()`](https://jsanchezalv.github.io/WARDEN/reference/rcond_gompertz_lu.md)
  : Draw from a Conditional Gompertz distribution (lower and upper
  bound)
- [`rdirichlet()`](https://jsanchezalv.github.io/WARDEN/reference/rdirichlet.md)
  : Draw from a dirichlet distribution based on number of counts in
  transition. Adapted from brms::rdirichlet
- [`rdirichlet_prob()`](https://jsanchezalv.github.io/WARDEN/reference/rdirichlet_prob.md)
  : Draw from a dirichlet distribution based on mean transition
  probabilities and standard errors
- [`pcond_gompertz()`](https://jsanchezalv.github.io/WARDEN/reference/pcond_gompertz.md)
  : Survival Probaility function for conditional Gompertz distribution
  (lower bound only)

## Other Functions and Utilities

- [`disc_cycle_v()`](https://jsanchezalv.github.io/WARDEN/reference/disc_cycle_v.md)
  : Cycle discounting for vectors
- [`disc_instant_v()`](https://jsanchezalv.github.io/WARDEN/reference/disc_instant_v.md)
  : Calculate instantaneous discounted costs or qalys for vectors
- [`disc_ongoing_v()`](https://jsanchezalv.github.io/WARDEN/reference/disc_ongoing_v.md)
  : Calculate discounted costs and qalys between events for vectors
- [`ast_as_list()`](https://jsanchezalv.github.io/WARDEN/reference/ast_as_list.md)
  : Transform a substituted expression to its Abstract Syntax Tree (AST)
  as a list
- [`extract_elements_from_list()`](https://jsanchezalv.github.io/WARDEN/reference/extract_elements_from_list.md)
  : Extracts items and events by looking into assignments, modify_event
  and new_event
- [`extract_from_reactions()`](https://jsanchezalv.github.io/WARDEN/reference/extract_from_reactions.md)
  : Extract all items and events and their interactions from the event
  reactions list
- [`extract_psa_result()`](https://jsanchezalv.github.io/WARDEN/reference/extract_psa_result.md)
  : Extract PSA results from a treatment
- [`tte.df`](https://jsanchezalv.github.io/WARDEN/reference/tte.df.md) :
  Example TTE IPD data

## Deprecated Functions

- [`add_item2()`](https://jsanchezalv.github.io/WARDEN/reference/add_item2.md)
  : Define parameters that may be used in model calculations (uses
  expressions)
- [`modify_item()`](https://jsanchezalv.github.io/WARDEN/reference/modify_item.md)
  : Modify the value of existing items
- [`modify_item_seq()`](https://jsanchezalv.github.io/WARDEN/reference/modify_item_seq.md)
  : Modify the value of existing items
- [`disc_cycle()`](https://jsanchezalv.github.io/WARDEN/reference/disc_cycle.md)
  : Cycle discounting
- [`disc_instant()`](https://jsanchezalv.github.io/WARDEN/reference/disc_instant.md)
  : Calculate instantaneous discounted costs or qalys
- [`disc_ongoing()`](https://jsanchezalv.github.io/WARDEN/reference/disc_ongoing.md)
  : Calculate discounted costs and qalys between events
