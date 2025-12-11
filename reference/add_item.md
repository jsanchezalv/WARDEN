# Define or append model inputs

Build a single [`{}`](https://rdrr.io/r/base/Paren.html) expression that
defines inputs for a simulation.

- Named args in `...` become assignments (`name <- expr`), e.g.,
  `add_item(a=5)`

- Unnamed args are inserted raw/unevaluated. If an unnamed arg is a
  [`{}`](https://rdrr.io/r/base/Paren.html) block, its statements are
  spliced (flattened). `add_item(pick_val_v(...))`

- Works with magrittr pipes: a leading `.` (the LHS) is resolved to its
  value; if that value is a [`{}`](https://rdrr.io/r/base/Paren.html)
  block (or list of expressions), it becomes the starting block.

- input argument can be used to handle alternative `add_item2` method,
  e.g. `add_item(input = {a <- 5})`

## Usage

``` r
add_item(..., .data = NULL, input)
```

## Arguments

- ...:

  Unevaluated arguments. Named → `name <- expr`; unnamed → raw expr.

- .data:

  Optional named argument: an existing
  [`{}`](https://rdrr.io/r/base/Paren.html) block (or list of
  expressions) to start from.

- input:

  Optional unevaluated expression or
  [`{}`](https://rdrr.io/r/base/Paren.html) block to splice in.

## Value

A single [`{}`](https://rdrr.io/r/base/Paren.html) call (language
object) ready for `load_inputs()`.

## Examples

``` r
library(magrittr)

add_item(input = {fl.idfs <-  0})
#> {
#>     fl.idfs <- 0
#> }
add_item(input = {
 util_idfs <- if(psa_bool){rnorm(1,0.8,0.2)} else{0.8}
 util.mbc <- 0.6
 cost_idfs <- 2500})
#> {
#>     util_idfs <- if (psa_bool) {
#>         rnorm(1, 0.8, 0.2)
#>     }
#>     else {
#>         0.8
#>     }
#>     util.mbc <- 0.6
#>     cost_idfs <- 2500
#> }
common_inputs <- add_item(input = {
pick_val_v(
  base      = l_statics[["base"]],
  psa       = pick_psa(
    l_statics[["function"]],
    l_statics[["n"]],
    l_statics[["a"]],
    l_statics[["b"]]
  ),
  sens      = l_statics[[sens_name_used]],
  psa_ind   = psa_bool,
  sens_ind  = sensitivity_bool,
  indicator = indicators_statics,
  names_out = l_statics[["parameter_name"]],
  deploy_env = TRUE #Note this option must be active if using it at add_item2
)
}
) 
```
