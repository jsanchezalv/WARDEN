### Function index by file

**`R/run_sim.R`** — Main user-facing simulation runner
- `run_sim()` [L150–523]: Run a DES (standard or constrained) across arms/sensitivities/PSA

**`R/run_sim_parallel.R`** — Parallel simulation runner
- `run_sim_parallel()` [L164–543]: Like `run_sim()` but uses `future`/`doFuture` for parallel execution

**`R/run_engine.R`** — Internal standard DES engine
- `run_engine()` [L20–290]: Internal engine processing events per patient per arm

**`R/run_engine_constrained.R`** — Internal constrained DES engine
- `run_engine_constrained()` [L47–386]: Internal engine for resource-constrained DES

**`R/engine_helper_f.R`** — Internal engine utilities
- `load_inputs()` [L41–44]: Evaluate and load unevaluated input expressions into environment
- `debug_inputs()` [L60–91]: Compare old vs new inputs for debug logging
- `initiate_evt()` [L109–130]: Initialise events for a patient at arm start
- `react_evt()` [L153–194]: Process the next event and apply reactions
- `eval_reactevt()` [L216–277]: Evaluate reaction list for a given event
- `get_input()` [L300–332]: Retrieve input value with type coercion
- `interval_out()` [L356–365]: Format simulation output as confidence interval string
- `.set_last_ctx()` [L367–393]: Internal — update error context beacon
- `log_add()` [L396–402]: Append an entry to the debug log
- `on_error_check()` [L404–416]: Wrap expression for error handling/continue-on-error
- `with_write_flags_lang()` [L467–512] *(exported)*: Track assignments in expression for debug mode
- `transform_debug()` [L527–542]: Reshape raw debug data into display format
- `export_log()` [L557–613]: Write debug log to file
- `expand_evts_bwd()` [L627–683]: Expand event data to fill time points (backward accumulation)
- `expand_evts_fwd()` [L697–753]: Expand event data to fill time points (forward accumulation)
- `compute_outputs_timseq()` [L785–996]: Compute timed-frequency outputs
- `compute_outputs()` [L1035–1392]: Compute all discounted outputs for a patient

**`R/input_f.R`** — Input construction, event queue wrappers, and utilities
- `replicate_profiles()` [L31–47]: Replicate patient profiles across arms
- `create_indicators()` [L67–94]: Create sensitivity/scenario indicator variables
- `sens_iterator()` [L115–120]: Iterate through DSA/scenario combinations
- `pick_psa()` [L192–203]: Draw PSA samples from a distribution function
- `pick_val_v()` [L265–357]: Select values across base/sensitivity/PSA with a vector of options
- `input_block()` [L439–531]: Build a `{}` block calling `pick_val_v`; `binary=TRUE` = 0/1 indicators per param, `binary=FALSE` (default) = grouped integer DSA labels via `dsa_indicators`
- `add_item()` [L534–590]: Add named items to the input list (does not support native pipe `|>`)
- `add_item2()` [L552–566] *(deprecated)*: Older item-adding variant (now merged into `add_item()`)
- `modify_item()` [L521–536] *(deprecated)*: Modify items in input list (no longer needed)
- `modify_item_seq()` [L572–597] *(deprecated)*: Sequential item modification (no longer needed)
- `queue_create()` [L610–612]: Create Rcpp-backed priority event queue
- `new_event()` [L639–675]: Add event(s) for a patient to the queue
- `next_event()` [L686–689]: Peek at the next event(s) in the queue
- `next_event_pt()` [L701–706]: Peek at next event(s) for a specific patient
- `pop_event()` [L716–719]: Remove the top event from the queue
- `pop_and_return_event()` [L1007–1009]: Remove and return the top event (creates new named list)
- `pop_into()` [L1012–1014]: Internal — pop into a pre-allocated 3-element list buffer (avoids list allocation)
- `remove_event()` [L743–779]: Remove named event(s) for a patient
- `modify_event()` [L806–848]: Modify event time/name for a patient (creates if missing by default)
- `queue_empty()` [L857–860]: Check if queue is empty
- `queue_size()` [L869–872]: Get number of events in queue
- `has_event()` [L883–887]: Check if a patient has a named event
- `get_event()` [L898–902]: Get time of a named event for a patient
- `resource_discrete()` [L960–1180]: Create Rcpp-backed discrete resource object; now accepts `discipline` ("FIFO"/"LIFO") and `max_queue`; new methods: `queue_wait_time`, `had_to_queue`, `time_in_use`, `utilization`, `n_using`, `total_patients_blocked`, `total_patients_queued`, `batch_seize`
- `print.resource_discrete()` [L1188–1196]: Print method for `resource_discrete`
- `seize()` [~L1200]: Acquire a resource for current patient; returns TRUE/FALSE/NA
- `release()` [~L1215]: Free a resource + auto-trigger next queued patient via `resume_event`
- `seize_all()` [~L1245]: Atomically seize multiple resources (C++); policies: all_or_none/sequential
- `release_all()` [~L1556]: Free multiple resources + purge from queues (C++) + schedule resume events only for freed resources
- `release_all_if_using()` [~L1598]: Free multiple resources only if using (no queue removal) + schedule resume events
- `shared_incr()` [~L1300]: Increment a `shared_input` counter, return new value
- `shared_decr()` [~L1310]: Decrement a `shared_input` counter, return new value
- `shared_input()` [~L1330–1400]: Create shared input object for constrained DES
- `add_reactevt()` [L1376–1399]: Add reactions (event-triggered logic) to the model
- `random_stream()` [L1439–1465]: Create a random number stream object for reproducibility
- `add_tte()` [L1494–1527]: Add time-to-event draws to the input list
- `adj_val()` [L1573–1613]: Integrate a value over a time interval (supports discounting)
- `disc_ongoing()` [L1635–1651]: Discount ongoing (flow) value between two time points
- `disc_instant()` [L1671–1681]: Discount an instantaneous value at a single time point
- `disc_cycle()` [L1709–1764]: Discount a cycle-based value
- `extract_from_reactions()` [L1799–1805]: Extract item/event references from reaction list
- `ast_as_list()` [L1861–1877]: Convert R AST expression to a nested list
- `extract_elements_from_list()` [L~2000+]: Extract assignments/references from AST list

**`R/calculator_f.R`** — Statistical/distributional utilities
- `draw_tte()` [L45–76]: Draw time-to-event from common parametric distributions (via flexsurv)
- `rdirichlet()` [L96–124]: Sample from a Dirichlet distribution
- `rdirichlet_prob()` [L148–176]: Sample Dirichlet with SE-based parameterisation
- `rbeta_mse()` [L197–208] / `qbeta_mse()` [L226–232]: Beta distribution parameterised by mean and SE
- `rgamma_mse()` [L253–267] / `qgamma_mse()` [L287–297]: Gamma distribution parameterised by mean and SE
- `rcond_gompertz_lu()` [L317–326] / `rcond_gompertz()` [L345–352] / `pcond_gompertz()` [L371–376]: Conditional Gompertz functions
- `rpoisgamma()` [L415–494]: Draw from Poisson-Gamma (negative binomial) distribution
- `cond_mvn()` [L522–564]: Conditional multivariate normal distribution
- `cond_dirichlet()` [L588–616]: Conditional Dirichlet distribution
- `discrete_resource_clone()` [L735–744]: Clone a `resource_discrete` object

**`R/RcppExports.R`** — Auto-generated Rcpp wrappers (do not edit manually)
- `luck_adj()` [L74–76]: Luck adjustment for correlated survival draws
- `qcond_exp()` [L91–93] / `qcond_weibull()` [L108–110] / `qcond_weibullPH()` [L125–127] / `qcond_llogis()` [L142–144] / `qcond_gompertz()` [L159–161] / `qcond_lnorm()` [L180–182] / `qcond_norm()` [L201–203] / `qcond_gamma()` [L222–224]: Vectorised conditional quantile functions
- `qtimecov()` [L378–380]: Time-to-event with time-varying covariates
- `disc_cycle_v()` [L422–424] / `disc_instant_v()` [L439–441] / `disc_ongoing_v()` [L457–459]: Vectorised Rcpp discounting functions
- Queue/resource C++ wrappers [L461–603]: `queue_create_cpp`, `new_event_cpp`, etc.

**`R/results_summary_f.R`** — Results summarisation
- `summary_results_det()` [L52–130]: Summarise deterministic simulation results
- `summary_results_sim()` [L160–242]: Summarise PSA simulation results
- `summary_results_sens()` [L273–371]: Summarise sensitivity/scenario analysis results
- `extract_psa_result()` [L401–406]: Extract a specific element from PSA results list
- `ceac_des()` [L445–482]: Compute cost-effectiveness acceptability curve (CEAC) data
- `evpi_des()` [L518–558]: Compute expected value of perfect information (EVPI)

**`R/old_R_nowincpp.R`** — Legacy R implementations (superseded by Rcpp, kept for reference)
- `luck_adj_old()` [L74–98] / `disc_ongoing_v_old()` [L117–130] / `disc_instant_v_old()` [L146–152] / `disc_cycle_v_old()` [L196–243]: Old vectorised R versions
- `qcond_gompertz_old()` [L260–268] / `qcond_exp_old()` [L285–291] / `qcond_weibull_old()` [L308–314] / `qcond_weibullPH_old()` [L330–340] / `qcond_llogis_old()` [L357–363] / `qcond_lnorm_old()` [L383–389] / `qcond_norm_old()` [L410–416] / `qcond_gamma_old()` [L437–443]: Old conditional quantile functions
- `qtimecov_old()` [L602–698]: Old time-varying covariate TTE function

**`R/data.R`** — Dataset documentation
- `tte.df`: Example time-to-event data frame for package examples