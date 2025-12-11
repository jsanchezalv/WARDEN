# Define parameters that may be used in model calculations (uses expressions)

Define parameters that may be used in model calculations (uses
expressions)

## Usage

``` r
add_item2(.data = NULL, input)
```

## Arguments

- .data:

  Existing data

- input:

  Items to define for the simulation as an expression (i.e., using )

## Value

A substituted expression to be evaluated by engine

## Details

DEPRECATED (old description): The functions to add/modify events/inputs
use named vectors or lists. If chaining together add_item2, it will just
append the expressions together in the order established.

If using `pick_val_v`, note it should be used with the
`deploy_env = TRUE` argument so that add_item2 process it correctly.
