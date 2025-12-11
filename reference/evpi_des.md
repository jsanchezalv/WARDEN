# Calculate the Expected Value of Perfect Information (EVPI) for a DES model with a PSA result

Calculate the Expected Value of Perfect Information (EVPI) for a DES
model with a PSA result

## Usage

``` r
evpi_des(wtp, results, interventions = NULL, sensitivity_used = 1)
```

## Arguments

- wtp:

  Vector of length \>=1 with the willingness to pay

- results:

  The list object returned by
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

- interventions:

  A character vector with the names of the interventions to be used for
  the analysis

- sensitivity_used:

  Integer signaling which sensitivity analysis to use

## Value

A data frame with the EVPI results

## Examples

``` r
res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
), merged_df = list(simulation = 1L, sensitivity = 1L))))

evpi_des(seq(from=10000,to=500000,by=10000),res)
#> # A tibble: 100 × 2
#>       wtp  evpi
#>     <dbl> <dbl>
#>  1  10000     0
#>  2  20000     0
#>  3  30000     0
#>  4  40000     0
#>  5  50000     0
#>  6  60000     0
#>  7  70000     0
#>  8  80000     0
#>  9  90000     0
#> 10 100000     0
#> # ℹ 90 more rows
```
