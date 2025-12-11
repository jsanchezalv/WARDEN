# Deterministic results for a specific treatment

Deterministic results for a specific treatment

## Usage

``` r
summary_results_det(out = results[[1]][[1]], arm = NULL, wtp = 50000)
```

## Arguments

- out:

  The final_output data frame from the list object returned by
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

- arm:

  The reference treatment for calculation of incremental outcomes

- wtp:

  Willingness to pay to have INMB

## Value

A dataframe with absolute costs, LYs, QALYs, and ICER and ICUR for each
intervention

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


summary_results_det(res[[1]][[1]],arm="int")
#>                        int     noint
#> costs             49921.64  41225.25
#> dcosts                0.00   8696.38
#> lys                   9.05      9.05
#> dlys                  0.00      0.00
#> qalys                 6.21      6.18
#> dqalys                0.00      0.03
#> ICER                    NA       Inf
#> ICUR                    NA 330825.35
#> INMB                    NA  -7382.03
#> costs_undisc      59831.36  49293.10
#> dcosts_undisc         0.00  10538.25
#> lys_undisc           10.90     10.90
#> dlys_undisc           0.00      0.00
#> qalys_undisc          7.50      7.47
#> dqalys_undisc         0.00      0.03
#> ICER_undisc             NA       Inf
#> ICUR_undisc             NA 389864.98
#> INMB_undisc             NA  -9186.73
#> c_default         49921.64  41225.25
#> dc_default            0.00   8696.38
#> c_default_undisc  59831.36  49293.10
#> dc_default_undisc     0.00  10538.25
#> q_default             6.21      6.18
#> dq_default            0.00      0.03
#> q_default_undisc      7.50      7.47
#> dq_default_undisc     0.00      0.03
```
