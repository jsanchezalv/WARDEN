# Example TTE IPD data

An example of TTE IPD data for the example_ipd file

## Usage

``` r
tte.df
```

## Format

### `tte.df`

A data frame with 1000 rows and 8 columns:

- USUBJID:

  Patient ID

- ARMCD, ARM:

  Arm code and variables

- PARAMCD, PARAM:

  Parameter

- AVAL, AVALCD:

  Values of interest

- CNSR:

  Censored observation?

## Source

Simulated through FlexsurvPlus package using sim_adtte(seed = 821, rho =
0, beta_1a = log(0.6), beta_1b = log(0.6), beta_pd = log(0.2))
