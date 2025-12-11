# Helper function to create a list with random draws or whenever a series of functions needs to be called. Can be implemented within `pick_val_v`.

Helper function to create a list with random draws or whenever a series
of functions needs to be called. Can be implemented within `pick_val_v`.

## Usage

``` r
pick_psa(f, ...)
```

## Arguments

- f:

  A string or vector of strings with the function to be called, e.g.,
  "rnorm"

- ...:

  parameters to be passed to the function (e.g., if "rnorm", arguments
  `n`, `mean`, `sd`)

## Value

List with length equal to `f` of parameters called

## Details

This function can be used to pick values for the PSA within
`pick_val_v.`

The function will ignore NA items within the respective parameter (see
example below). If an element in f is NA (e.g., a non PSA input) then it
will return NA as its value This feature is convenient when mixing
distributions with different number of arguments, e.g., `rnorm` and
`rgengamma`.

While it's slightly lower than individually calling each function, it
makes the code easier to read and more transparent

## Examples

``` r
params <- list(
param=list("a","b"),
dist=list("rlnorm","rnorm"),
n=list(4,1),
a=list(c(1,2,3,4),1),
b=list(c(0.5,0.5,0.5,0.5),0.5),
dsa_min=list(c(1,2,3,4),2),
dsa_max=list(c(1,2,3,4),3)
)
pick_psa(params[["dist"]],params[["n"]],params[["a"]],params[["b"]])
#> [[1]]
#> [1]  1.675586 10.981649 49.775128 34.175660
#> 
#> [[2]]
#> [1] 1.43065
#> 

#It works with functions that require different number of parameters
params <- list(
 param=list("a","b","c"),
 dist=list("rlnorm","rnorm","rgengamma"),
 n=list(4,1,1),
 a=list(c(1,2,3,4),1,0),
 b=list(c(0.5,0.5,0.5,0.5),0.5,1),
 c=list(NA,NA,0.2),
 dsa_min=list(c(1,2,3,4),2,1),
 dsa_max=list(c(1,2,3,4),3,3)
)

pick_psa(params[["dist"]],params[["n"]],params[["a"]],params[["b"]],params[["c"]])
#> [[1]]
#> [1]  3.709214  2.529025 25.594074 43.656966
#> 
#> [[2]]
#> [1] 1.273049
#> 
#> [[3]]
#> [1] 1.151862
#> 

#Can be combined with multiple type of functions and distributions if parameters are well located

params <- list(
param=list("a","b","c","d"),
dist=list("rlnorm","rnorm","rgengamma","draw_tte"),
n=list(4,1,1,1),
a=list(c(1,2,3,4),1,0,"norm"),
b=list(c(0.5,0.5,0.5,0.5),0.5,1,1),
c=list(NA,NA,0.2,0.5),
c=list(NA,NA,NA,NA), #NA arguments will be ignored
dsa_min=list(c(1,2,3,4),2,1,0),
dsa_max=list(c(1,2,3,4),3,3,2)
)
```
