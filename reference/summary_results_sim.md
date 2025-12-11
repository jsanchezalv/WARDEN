# Summary of PSA outputs for a treatment

Summary of PSA outputs for a treatment

## Usage

``` r
summary_results_sim(out = results[[1]], arm = NULL, wtp = 50000)
```

## Arguments

- out:

  The output_sim data frame from the list object returned by
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

- arm:

  The reference treatment for calculation of incremental outcomes

- wtp:

  Willingness to pay to have INMB

## Value

A data frame with mean and 95% CI of absolute costs, LYs, QALYs, ICER
and ICUR for each intervention from the PSA samples

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


summary_results_sim(res[[1]],arm="int")
#>                                                int
#> costs                      49,922 (49,922; 49,922)
#> dcosts                                    0 (0; 0)
#> lys                              9.05 (9.05; 9.05)
#> dlys                                      0 (0; 0)
#> qalys                            6.21 (6.21; 6.21)
#> dqalys                                    0 (0; 0)
#> ICER                                  NaN (NA; NA)
#> ICUR                                  NaN (NA; NA)
#> INMB                                  NaN (NA; NA)
#> costs_undisc               59,831 (59,831; 59,831)
#> dcosts_undisc                             0 (0; 0)
#> lys_undisc                       10.9 (10.9; 10.9)
#> dlys_undisc                               0 (0; 0)
#> qalys_undisc                        7.5 (7.5; 7.5)
#> dqalys_undisc                             0 (0; 0)
#> ICER_undisc                           NaN (NA; NA)
#> ICUR_undisc                           NaN (NA; NA)
#> INMB_undisc                           NaN (NA; NA)
#> c_default         49,921.64 (49,921.64; 49,921.64)
#> dc_default                                0 (0; 0)
#> c_default_undisc  59,831.36 (59,831.36; 59,831.36)
#> dc_default_undisc                         0 (0; 0)
#> q_default                        6.21 (6.21; 6.21)
#> dq_default                                0 (0; 0)
#> q_default_undisc                    7.5 (7.5; 7.5)
#> dq_default_undisc                         0 (0; 0)
#>                                              noint
#> costs                      41,225 (41,225; 41,225)
#> dcosts                        8,696 (8,696; 8,696)
#> lys                              9.05 (9.05; 9.05)
#> dlys                                      0 (0; 0)
#> qalys                            6.18 (6.18; 6.18)
#> dqalys                        0.026 (0.026; 0.026)
#> ICER                                Inf (Inf; Inf)
#> ICUR                    330,825 (330,825; 330,825)
#> INMB                       -7,382 (-7,382; -7,382)
#> costs_undisc               49,293 (49,293; 49,293)
#> dcosts_undisc              10,538 (10,538; 10,538)
#> lys_undisc                       10.9 (10.9; 10.9)
#> dlys_undisc                               0 (0; 0)
#> qalys_undisc                     7.47 (7.47; 7.47)
#> dqalys_undisc                 0.027 (0.027; 0.027)
#> ICER_undisc                         Inf (Inf; Inf)
#> ICUR_undisc             389,865 (389,865; 389,865)
#> INMB_undisc                -9,187 (-9,187; -9,187)
#> c_default         41,225.25 (41,225.25; 41,225.25)
#> dc_default        8,696.381 (8,696.381; 8,696.381)
#> c_default_undisc     49,293.1 (49,293.1; 49,293.1)
#> dc_default_undisc 10,538.25 (10,538.25; 10,538.25)
#> q_default                        6.18 (6.18; 6.18)
#> dq_default                    0.026 (0.026; 0.026)
#> q_default_undisc                 7.47 (7.47; 7.47)
#> dq_default_undisc             0.027 (0.027; 0.027)
```
