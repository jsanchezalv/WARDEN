# WARDEN 1.3.4
* Unit tests for shared_inputs and run_sim_parallel added.
* add_item and add_item2 have now been integrated into add_item, with both behaviors being accepted. This means load_inputs has been overwritten with the load_inputs2 form.

# WARDEN 1.3.3
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
