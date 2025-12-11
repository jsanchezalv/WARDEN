# Modify the value of existing items

Modify the value of existing items

## Usage

``` r
modify_item(list_item)
```

## Arguments

- list_item:

  A list of items and their values or expressions

## Value

No return value, modifies/adds item to the environment and integrates it
with the main list for storage

## Details

DEPRECATED (old description): The functions to add/modify events/inputs
use lists. Whenever several inputs/events are added or modified, it's
recommended to group them within one function, as it reduces the
computation cost. So rather than use two `modify_item` with a list of
one element, it's better to group them into a single `modify_item` with
a list of two elements.

Note that `modify_item` nor `modify_item_seq` can work on subelements
(e.g., `modify_item(list(obj$item = 5))` will not work as intended, for
that is better to assign directly using the expression approach, so
`obj$item <- 5`).

Costs and utilities can be modified by using the construction
`type_name_category`, where type is either "qaly" or "cost", name is the
name (e.g., "default") and category is the category used (e.g.,
"instant"), so one could pass `cost_default_instant` and modify the
cost. This will overwrite the value defined in the corresponding
cost/utility section.

This function is intended to be used only within the `add_reactevt`
function in its `input` parameter and should not be run elsewhere or it
will return an error.

## Examples

``` r
add_reactevt(name_evt = "idfs",input = {modify_item(list("cost.it"=5))})
#> $idfs
#> $idfs$react
#> {
#>     modify_item(list(cost.it = 5))
#> }
#> 
#> 
```
