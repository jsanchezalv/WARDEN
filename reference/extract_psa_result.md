# Extract PSA results from a treatment

Extract PSA results from a treatment

## Usage

``` r
extract_psa_result(x, element)
```

## Arguments

- x:

  The output_sim data frame from the list object returned by
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

- element:

  Variable for which PSA results are being extracted (single string)

## Value

A dataframe with PSA results from the specified intervention

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


extract_psa_result(res[[1]],"total_costs")
#>        int    noint simulation     element
#> 1 49921.64 41225.25          1 total_costs
```
