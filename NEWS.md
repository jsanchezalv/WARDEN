# RDICE 0.96
* Update based on validation comments from Gabriel. 
* Renamed conditional quantile functions for consistency.

# RDICE 0.95
* Seeds used by default have been changed to guarantee uniqueness
* Added possibility of continuing to the next simulation on error (if it occurs at the patient/arm level, not at the statics/structural loading level)
* Debug mode now exports log even if the simulation stops due to error. If combined with continue on error, it will continue
to export log with the timestamp

# RDICE 0.94
* Added possibility of accumulating outputs continuously backward or forward using the accum_backward option in run_sim and run_sim_parallel

# RDICE 0.93
* Conditional quantile functions added and adjusted
* Luck adjustment function added with instructions for user

# RDICE 0.92
* Progress bar added for both parallel and standard computing model
* To use progress bar while in batch mode for a quarto document, make sure to add in the knitr options 
knitr:
  opts_chunk:
    R.options:
      progressr.enable: true

# RDICE 0.91
* Debug mode exports a txt file

# RDICE 0.9
* Warning, this commit will change previous results. 
* Sensitivity-level and simulaton-level seeds have been moved outside the input loading loop as this caused correlation between inputs loaded at those stages. 

# RDICE 0.5

* Initial set-up of news file
* Summary of inputs has been overhauled to provide INMB through a WTP argument
* Summary now can also be provided across analyses to quickly obtain DSA/scenario analysis results summarized
