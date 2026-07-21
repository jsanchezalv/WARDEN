# WARDEN 2.0.6

* `add_reactevt()` now detects when a character string is accidentally passed as `.data` (the pipe-first argument) and errors with guidance to use `name_evt =` explicitly.
* `draw_categorical()` and `qcategorical()` are new utility functions that return a 1-based integer index from a vector of probabilities, replacing the error-prone `findInterval(luck, cumsum(probs)) + 1L` pattern.
* `print.warden_results()` is a new print method for results objects returned by `run_sim()` and `run_sim_parallel()`, showing the number of analyses, simulations, arms, and access-pattern hints.
* `random_stream()` gains a `strict` argument (default `FALSE`). When `TRUE`, `draw_n()` throws an error on stream exhaustion instead of auto-regenerating with a warning.
* `run_sim()` and `run_sim_parallel()` now check for integer overflow in seed computation at call time, before the simulation starts.
* `run_sim()` and `run_sim_parallel()` now validate that outcome accumulator variables (e.g., `q_default`, `c_default`) are scalar (length 1, unnamed) at the start of the first patient-arm evaluation.
* `scale_remaining_time()` is a new utility that wraps the common `curtime + (get_event(x) - curtime) / hr` pattern for rescaling event times by a hazard ratio.
* `seize()` gains a `priority` argument (default `1L`) forwarded to `attempt_block()`. Passing `seize_all()`-specific arguments (e.g., `accum_queue`, `force_unblock`) now produces a clear error pointing to `seize_all()`.
* `summary_results_det()`, `summary_results_sim()`, and `summary_results_sens()` now auto-detect the nesting level of the results object and either auto-unwrap with a message or error with guidance.
* `validate_model()` is a new function that performs static checks on model inputs (arm coverage, event/reaction matching, seed overflow, dimensions) before running the simulation.
* `waning_hr()` is a new utility function that computes an effective hazard ratio under linear, exponential, or instant treatment waning, implementing a common NICE submission pattern.

# WARDEN 2.0.5
* `release_all()` and `release_all_if_using()` now deduplicate resume events: when the same patient is first in queue on multiple released resources, only one resume event is scheduled (prevents double-seize on retry).
* `seize_all()` gains an `accum_queue` argument (default `TRUE`). When `FALSE`, retries do not accumulate phantom queue entries on bottleneck resources — existing entries serve as placeholders.
* `seize_all()` no longer rejects (returns `NA`) when retrying on resources with `allow_multiple_queue = FALSE`. Previously, a retry was incorrectly rejected because the patient was already queued; now the existing queue entry is recognized and preserved.
* `new_event()` and `modify_event()` now reject `NA` and `NaN` event times with an informative error; `Inf` event times remain valid.
* `release()` and `attempt_free()` now enforce indivisible units: `release(amount = NULL)` (the default) releases one unit (equivalent to `amount = 1`), and `release(amount = X)` requires an exact match with the seized amount — mismatches throw an error. `remove_all = TRUE` still removes all usage entries regardless of amount.
* `release_all()` and `release_all_if_using()` now follow an all-or-none policy: both functions only act if the patient is currently using **all** listed resources. If the patient holds some but not all, they do nothing and emit an informative message. `release_all(amounts = NULL)` purges queue entries as before; `release_all(amounts = c(...))` performs indivisible per-resource release without touching queue entries. `release_all()` now also purges queue entries when the patient is not using any of the listed resources (e.g. a patient who queued via `seize_all()` and then died before acquiring); previously, calling `release_all()` at a patient's death event left dead patients stranded in resource queues, causing them to be rescheduled when another patient later released.
* `release_if_using()` is a new exported function that releases a resource for the current patient only if they are currently using it, doing nothing otherwise. Unlike `release()`, it never errors on a non-using patient and never removes queue entries.
* `resource_discrete()` gains an `allow_multiple_queue` argument (default `TRUE`). When `FALSE`, a patient that already has a queued entry will be rejected (`attempt_block()` returns `NA`) on a second seize attempt rather than being allowed to queue again.
* `seize_all()` with `policy = "all_or_none"` now queues the patient on **all** bottleneck resources simultaneously (previously only the first unavailable resource was queued). If any bottleneck would reject queuing, no resource is queued and `NA` is returned. `seize_all()` also gains a `force_unblock` argument: when `TRUE` and all resources have capacity but the patient is blocked by queue position on some resources, the patient is forcibly moved to the front of those queues and acquires all resources.

# WARDEN 2.0.4
* `release()` is a new helper that frees a `resource_discrete` for the current patient and, when `resume_event` is supplied, automatically schedules that event for the next patient in the resource queue — eliminating the manual `if(success & queue_size > 0) new_event(...)` pattern.
* `release_all()` is a new helper that fully purges the current patient from a list of `resource_discrete` objects (removes from using and from all queue entries) and optionally schedules per-resource resume events. Resume events are only triggered for resources where the patient was actually using (capacity freed).
* `release_all_if_using()` is a new helper that frees the current patient from a list of `resource_discrete` objects only if they are currently using them. Does not touch queue entries, preserving queue position for multi-resource workflows.
* `resource_discrete()` gains `discipline` (`"FIFO"` or `"LIFO"`) and `max_queue` arguments. `discipline` controls ordering within the same priority level. `max_queue` caps the waiting list so that `attempt_block()` returns `NA` (patient rejected) when the queue is full, enabling M/M/c/k systems.
* `resource_discrete()` gains new statistics methods: `queue_wait_time()`, `queue_wait_time_current()`, `had_to_queue()`, `time_in_use()`, `utilization()`, `n_using()`, `total_patients_blocked()`, and `total_patients_queued()`. `queue_wait_time_current()` returns `0` at queue entry, grows with elapsed time at each subsequent event while waiting, and returns the total wait on acquisition — replacing the need to manually track `time_start_queue`. `queue_wait_time()` returns the final stored wait only after dequeue.
* `resource_discrete()` gains `batch_seize()` for seizing a resource for multiple patients in a single C++ call, and `attempt_block()` / `attempt_free()` now accept an `amount` argument for multi-unit acquisitions by a single patient.
* `seize()` is a new thin wrapper around `resource$attempt_block()` that reads `i` and `curtime` from the calling environment, returning `TRUE` (acquired), `FALSE` (queued), or `NA` (rejected).
* `seize_all()` is a new helper that atomically seizes a list of resources in C++ under `"all_or_none"` or `"sequential"` policy, enabling deadlock-free surgery-scheduling and philosophers-problem patterns. Fixed a bug where `"all_or_none"` could partially acquire resources when one resource had free capacity but patients already in its queue.
* `shared_decr()` and `shared_incr()` are new one-liner helpers for `shared_input` counters that increment or decrement by `delta` and return the new value.

# WARDEN 2.0.3
* `add_item()` now works correctly with the native pipe (`|>`). The `.data` argument has been moved to the first position (`.data = NULL, ..., input`), so the LHS of `|>` is naturally routed to `.data` without relying on magrittr's `.` symbol. Existing code using `%>%`, `input=`, or named `...` arguments is unaffected (#TODO).
* `input_block()` and `run_sim()` now correctly handle multiple blocks spanning different simulation levels (e.g., `common_all_inputs` and `common_pt_inputs`): `n_sensitivity` is summed across all blocks, and binary-mode parameter offsets are injected automatically so each block activates at the right DSA iteration.
* `input_block()` is a new helper that builds a complete `pick_val_v()` expression from explicit `base`, `psa`, `sens`, and `names_out` arguments. The `binary` and `dsa_indicators` parameters have been renamed to `indicator_sens_binary` and `sens_indicators` respectively to align with `pick_val_v()`. The `dsa_names` argument now defaults to `NULL`, meaning all `sensitivity_names` are treated as scenarios (one iteration per name); supply `dsa_names` explicitly to designate which names are DSA directions. Setting an entry in `sens_indicators` to `0` now permanently excludes that parameter from variation in both DSA and scenario analyses, and the engine automatically deduces `n_sensitivity` from the number of active parameters/groups.
* `pick_val_v()` now correctly respects `indicator_psa` in grouped mode (`indicator_sens_binary = FALSE`). Previously, when `sens_bool = TRUE` and `psa_bool = TRUE`, all parameters drew from PSA regardless of `indicator_psa`; parameters with `indicator_psa = 0` now correctly draw from `base`.

# WARDEN 2.0.1
*adj_val now accepts a vectorized_f argument to speed computations in the case of vectorized functions

# WARDEN 2.0.0
* Rcpp event based handler has been created using queues for high efficiency. In the new system,
a unique event per patient is accepted in the queue, so small changes may be needed in codes where
multiple equally named events were set to accommodate the new system. This is to ensure best practices and clarity
regarding which events are modified, in the queue, etc.
* Simulation and engines have been redesigned and simplified to more cleanly handle and report errors, debug mode and continue on error functions.
* Constrained based engine has been created, which allows to run resource constrained DES, where resource
and/or inputs can be shared across patients within an arm. This can be activated by using run_sim with constrained = TRUE.
* Rcpp based discrete resources with their own queuing system have been implemented (only interacted through R,
working similar to R6 objects). This can be created by using resource_discrete()
* Shared inputs have been created to be used in constrained DES (works similar to R6 objects). This can be created by using shared_input().
* New vignette showcasing an example for constrained DES has been created.

# WARDEN 1.3.5
* Accumulation backwards now uses active bindings to recognize which ongoing outputs are being modified.
* Deprecated modify_item_seq, modify_item as they are no longer needed to run even if backwards = TRUE. Removed from all examples.

# WARDEN 1.3.4
* Unit tests for shared_inputs and run_sim_parallel added.
* add_item and add_item2 have now been integrated into add_item, with both behaviors being accepted. This means load_inputs has been overwritten with the load_inputs2 form.

# WARDEN 1.3.0
* Rcpp versions of key functions have been implemented for speed improvements:
disc_cycle_v, disc_ongoing_v, disc_instant_v, qcond_*, qtimecov, luck_adj

# WARDEN 1.2.4
* Added "adj_val" function, to reduce the need to call cycle events which can hinder model speed.
* Added "qtimecov" function, to predict time to events with time changing covariates, 
to reduce the need to call cycle events which can hinder model speed.

# WARDEN 1.2.3
* Added model run unit tests
* Now the model can be run without specifying utilities, costs or other outputs (only LYs will be accounted for)
* Fixed an issue when timed_freq was active (now outputs are correct)

# WARDEN 1.2.2
* Per cycle outcomes now correctly work, and a maximum number of cycles argument has been implemented.

# WARDEN 1.2.1
* Discounting now correctly allocates to drq or drc depending on output. Other inputs use drc.

# WARDEN 1.2.0
* Speed gains from rewritten internal compute_outputs function.
* Implemented "random_stream" function (using method similar to R6) to facilitate careful handling of random numbers.

# WARDEN 1.1.0
* Added "add_item2", which allows to incorporate expressions directly instead of a list, being faster and more consistent with "add_tte" and "add_react_evt".
* Engine now fully utilizes environment for slightly faster analysis when loading inputs.

# WARDEN 1.0
* Major update: now the engine uses environments instead of lists, which allows user to remove "modify_item" and "modify_item_seq" from their code,
improving running speed by 20-40%
* Secondary changes to accommodate this update applied throughout (extract reactions, debug mode).
* Debug mode now uses abstract syntax tree to capture assignments, which can be limited in the presence of dynamic code assignments.
* sens_iterator function added to facilitate looping through DSA or multiple scenarios in a single model run. 
See the input selectors vignette in the website for an example.

# WARDEN 0.99.3
* Minor fix in run_sim_parallel to ensure compatibility with future package (had some "=T" instead of "=TRUE")
* Added BugsReport link in Description

# WARDEN 0.99.2
* Added qgamma_mse function
* Added two articles for the website explaining more in detail WARDEN and the use of sobol sequences

# WARDEN 0.99.1
* CRAN feedback implemented, including changes in documentation:
* runif_stream function has been removed due to violation of CRAN policy of global environment modification.
It is suggested that the user employs different methods (e.g., pre-drawing random numbers)
* Now the user is allowed to select the starting seed for the analysis

# WARDEN 0.99
* Now debug and continue_on_error will work on all stages, not only for simulations
* CRAN preparation changes

# WARDEN 0.98
* Repository is now public, Github Website has been set up
* Website references are now split by topic
* Added auxiliary functions to extract items and events from reactions to easily see interconections in models

# WARDEN 0.97
* Update based on validation comments from Gabriel
* Modified conditional quantile functions weibull and llogistic to better match default R stats behavior
* Set License to be GPL >=3

# WARDEN 0.96
* Update based on validation comments from Gabriel. 
* Renamed conditional quantile functions for consistency.

# WARDEN 0.95
* Seeds used by default have been changed to guarantee uniqueness
* Added possibility of continuing to the next simulation on error (if it occurs at the patient/arm level, not at the statics/structural loading level)
* Debug mode now exports log even if the simulation stops due to error. If combined with continue on error, it will continue
to export log with the timestamp

# WARDEN 0.94
* Added possibility of accumulating outputs continuously backward or forward using the accum_backward option in run_sim and run_sim_parallel

# WARDEN 0.93
* Conditional quantile functions added and adjusted
* Luck adjustment function added with instructions for user

# WARDEN 0.92
* Progress bar added for both parallel and standard computing model
* To use progress bar while in batch mode for a quarto document, make sure to add in the knitr options 
knitr:
  opts_chunk:
    R.options:
      progressr.enable: true

# WARDEN 0.91
* Debug mode exports a txt file

# WARDEN 0.9
* Warning, this commit will change previous results. 
* Sensitivity-level and simulaton-level seeds have been moved outside the input loading loop as this caused correlation between inputs loaded at those stages. 

# WARDEN 0.5

* Initial set-up of news file
* Summary of inputs has been overhauled to provide INMB through a WTP argument
* Summary now can also be provided across analyses to quickly obtain DSA/scenario analysis results summarized
