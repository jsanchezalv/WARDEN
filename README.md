
# RDICE: Discrete Event Simulation for Cost-Effectiveness Modeling

## Introduction

`RDICE` is a user-friendly package that facilitates the use of discrete
event simulations without resource constraints for cost-effectiveness
analysis. The package supports a flexible, practical approach to
discrete event simulation while keeping an acceptable performance
through the use of parallel computing.

The current version supports:

- Discrete event simulation models, Markov/semi-Markov models and hybrid
  models using parallel and non-parallel engines
- Seamlessly integrating data.frames and other objects into the model
- Delayed execution of the main inputs to facilitate readability of the
  model
- Implementation of structural and parameter uncertainty
- Helper functions to facilitate drawing of time to events and the use
  of hazard ratios
- Performing cost-effectiveness and uncertainty analysis

It is recommended that the user checks the vignettes, first the simple
[Sick-Sicker-Dead
model](https://github.com/javiersanchezalvarezevidera/RDICE/articles/example_ssd.html)
and then the more complex model for [early breast
cancer](https://github.com/javiersanchezalvarezevidera/RDICE/articles/example_eBC.html).
The
[markov](https://github.com/javiersanchezalvarezevidera/RDICE/articles/example_markov.html)
example shows how to run a cohort Markov model while using the same
modeling framework. Similarly, a simulation based Markov model could be
run. Structural and parametric uncertainty are explored in the
[corresponding
vignette](https://github.com/javiersanchezalvarezevidera/RDICE/articles/example_uncertainty.html).
The [IPD
vignette](https://github.com/javiersanchezalvarezevidera/RDICE/articles/example_ipd.html)
shows how RDICE can be used when individual patient data is available.

## Documentation

Have a look at the [package home
site](https://github.com/javiersanchezalvarezevidera/RDICE/index.html)
for more details on documentation and specific tutorials.

For more details on the code, check our [Github
repository](https://github.com/javiersanchezalvarezevidera/RDICE).

## Installation

`RDICE` can the be installed directly from this repo via

``` r
# install.packages("devtools")
devtools::install_github("javiersanchezalvarezevidera/RDICE", ref="main")
```

## Citation

If you use `RDICE`, please contact the authors for the most up to date
appropiate citation.
