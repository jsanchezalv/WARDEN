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
