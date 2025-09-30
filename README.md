
# WARDEN: Workflows for health technology Assessments in R using Discrete EveNts

## Introduction

`WARDEN` is a user-friendly package that facilitates the use of discrete
event simulations with or without resource constraints for cost-effectiveness
analysis. The package supports a flexible, transparent, practical
approach to discrete event simulation while keeping an acceptable
performance. WARDEN now integrates R and C++ code to allow resource constrained
DES to be used with valid levels of performance.

The current version supports:

- Discrete event simulation models, Markov/semi-Markov models and hybrid
  models using parallel and non-parallel engines
- Seamlessly integrating data.frames and other objects into the model
- Delayed execution of the main inputs to facilitate readability of the
  model
- Implementation of structural and parameter uncertainty
- Helper functions to facilitate drawing of time to events and the use
  of hazard ratios, as well as other functions to facilitate
  transparency
- Performing cost-effectiveness and uncertainty analysis
- Resource constrained and unconstrained DES

It is recommended that the user checks the vignettes, first the simple
[Sick-Sicker-Dead
model](https://jsanchezalv.github.io/WARDEN/articles/example_ssd.html)
and then the more complex model for [early breast
cancer](https://jsanchezalv.github.io/WARDEN/articles/example_eBC.html).
The
[markov](https://jsanchezalv.github.io/WARDEN/articles/example_markov.html)
example shows how to run a cohort Markov model while using the same
modeling framework. Similarly, a simulation based Markov model could be
run. Structural and parametric uncertainty are explored in the
[corresponding
vignette](https://jsanchezalv.github.io/WARDEN/articles/example_uncertainty.html).
The [IPD
vignette](https://jsanchezalv.github.io/WARDEN/articles/example_ipd.html)
shows how WARDEN can be used when individual patient data is available. For
resource constrained, it's recommended to read its 
[vignette](https://jsanchezalv.github.io/WARDEN/articles/example_ssd_constrained.html).
There are some additional vignettes showcasing specific functions or functionalities.

## Documentation

Have a look at the [package home
site](https://jsanchezalv.github.io/WARDEN/index.html) for more details
on documentation and specific tutorials.

For more details on the code, check our [Github
repository](https://github.com/jsanchezalv/WARDEN), or our [CRAN
site](https://cran.r-project.org/package=WARDEN).

## Installation

`WARDEN` can now be installed directly from CRAN or this repo via

``` r
install.packages("WARDEN") #CRAN version 

# install.packages("devtools")
devtools::install_github("jsanchezalv/WARDEN", ref="main") #github version
```

## Citation

If you use `WARDEN`, please contact the authors for the most up to date
appropriate citation.

<!-- badges: start -->

[![R-CMD-check](https://github.com/jsanchezalv/WARDEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jsanchezalv/WARDEN/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# WARDEN <a href="https://jsanchezalv.github.io/WARDEN/"><img src="man/figures/logo.png" align="right" height="137" alt="WARDEN website" /></a>
