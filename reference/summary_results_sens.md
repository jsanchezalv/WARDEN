# Summary of sensitivity outputs for a treatment

Summary of sensitivity outputs for a treatment

## Usage

``` r
summary_results_sens(out = results, arm = NULL, wtp = 50000)
```

## Arguments

- out:

  The list object returned by
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

- arm:

  The reference treatment for calculation of incremental outcomes

- wtp:

  Willingness to pay to have INMB

## Value

A data frame with each sensitivity output per arm

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


summary_results_sens(res,arm="int")
#>        arm analysis analysis_name          variable
#>     <char>    <int>        <char>            <fctr>
#>  1:    int        1                           costs
#>  2:  noint        1                           costs
#>  3:    int        1                          dcosts
#>  4:  noint        1                          dcosts
#>  5:    int        1                             lys
#>  6:  noint        1                             lys
#>  7:    int        1                            dlys
#>  8:  noint        1                            dlys
#>  9:    int        1                           qalys
#> 10:  noint        1                           qalys
#> 11:    int        1                          dqalys
#> 12:  noint        1                          dqalys
#> 13:    int        1                            ICER
#> 14:  noint        1                            ICER
#> 15:    int        1                            ICUR
#> 16:  noint        1                            ICUR
#> 17:    int        1                            INMB
#> 18:  noint        1                            INMB
#> 19:    int        1                    costs_undisc
#> 20:  noint        1                    costs_undisc
#> 21:    int        1                   dcosts_undisc
#> 22:  noint        1                   dcosts_undisc
#> 23:    int        1                      lys_undisc
#> 24:  noint        1                      lys_undisc
#> 25:    int        1                     dlys_undisc
#> 26:  noint        1                     dlys_undisc
#> 27:    int        1                    qalys_undisc
#> 28:  noint        1                    qalys_undisc
#> 29:    int        1                   dqalys_undisc
#> 30:  noint        1                   dqalys_undisc
#> 31:    int        1                     ICER_undisc
#> 32:  noint        1                     ICER_undisc
#> 33:    int        1                     ICUR_undisc
#> 34:  noint        1                     ICUR_undisc
#> 35:    int        1                     INMB_undisc
#> 36:  noint        1                     INMB_undisc
#> 37:    int        1                       c_default
#> 38:  noint        1                       c_default
#> 39:    int        1                      dc_default
#> 40:  noint        1                      dc_default
#> 41:    int        1                c_default_undisc
#> 42:  noint        1                c_default_undisc
#> 43:    int        1               dc_default_undisc
#> 44:  noint        1               dc_default_undisc
#> 45:    int        1                       q_default
#> 46:  noint        1                       q_default
#> 47:    int        1                      dq_default
#> 48:  noint        1                      dq_default
#> 49:    int        1                q_default_undisc
#> 50:  noint        1                q_default_undisc
#> 51:    int        1               dq_default_undisc
#> 52:  noint        1               dq_default_undisc
#>        arm analysis analysis_name          variable
#>     <char>    <int>        <char>            <fctr>
#>                                value
#>                               <char>
#>  1:          49,922 (49,922; 49,922)
#>  2:          41,225 (41,225; 41,225)
#>  3:                         0 (0; 0)
#>  4:             8,696 (8,696; 8,696)
#>  5:                9.05 (9.05; 9.05)
#>  6:                9.05 (9.05; 9.05)
#>  7:                         0 (0; 0)
#>  8:                         0 (0; 0)
#>  9:                6.21 (6.21; 6.21)
#> 10:                6.18 (6.18; 6.18)
#> 11:                         0 (0; 0)
#> 12:             0.026 (0.026; 0.026)
#> 13:                     NaN (NA; NA)
#> 14:                   Inf (Inf; Inf)
#> 15:                     NaN (NA; NA)
#> 16:       330,825 (330,825; 330,825)
#> 17:                     NaN (NA; NA)
#> 18:          -7,382 (-7,382; -7,382)
#> 19:          59,831 (59,831; 59,831)
#> 20:          49,293 (49,293; 49,293)
#> 21:                         0 (0; 0)
#> 22:          10,538 (10,538; 10,538)
#> 23:                10.9 (10.9; 10.9)
#> 24:                10.9 (10.9; 10.9)
#> 25:                         0 (0; 0)
#> 26:                         0 (0; 0)
#> 27:                   7.5 (7.5; 7.5)
#> 28:                7.47 (7.47; 7.47)
#> 29:                         0 (0; 0)
#> 30:             0.027 (0.027; 0.027)
#> 31:                     NaN (NA; NA)
#> 32:                   Inf (Inf; Inf)
#> 33:                     NaN (NA; NA)
#> 34:       389,865 (389,865; 389,865)
#> 35:                     NaN (NA; NA)
#> 36:          -9,187 (-9,187; -9,187)
#> 37: 49,921.64 (49,921.64; 49,921.64)
#> 38: 41,225.25 (41,225.25; 41,225.25)
#> 39:                         0 (0; 0)
#> 40: 8,696.381 (8,696.381; 8,696.381)
#> 41: 59,831.36 (59,831.36; 59,831.36)
#> 42:    49,293.1 (49,293.1; 49,293.1)
#> 43:                         0 (0; 0)
#> 44: 10,538.25 (10,538.25; 10,538.25)
#> 45:                6.21 (6.21; 6.21)
#> 46:                6.18 (6.18; 6.18)
#> 47:                         0 (0; 0)
#> 48:             0.026 (0.026; 0.026)
#> 49:                   7.5 (7.5; 7.5)
#> 50:                7.47 (7.47; 7.47)
#> 51:                         0 (0; 0)
#> 52:             0.027 (0.027; 0.027)
#>                                value
#>                               <char>
```
