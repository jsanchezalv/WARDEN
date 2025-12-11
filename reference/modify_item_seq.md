# Modify the value of existing items

Modify the value of existing items

## Usage

``` r
modify_item_seq(...)
```

## Arguments

- ...:

  A list of items and their values or expressions. Will be evaluated
  sequentially (so one could have list(a= 1, b = a +2 ))

## Value

No return value, modifies/adds items sequentially and deploys to the
environment and with the main list for storage

## Details

The functions to add/modify events/inputs use lists. Whenever several
inputs/events are added or modified, it's recommended to group them
within one function, as it reduces the computation cost. So rather than
use two `modify_item` with a list of one element, it's better to group
them into a single `modify_item` with a list of two elements.

Note that `modify_item` nor `modify_item_seq` can work on subelements
(e.g., `modify_item_seq(list(obj$item = 5))` will not work as intended,
for that is better to assign directly using the expression approach, so
`obj$item <- 5`).

Costs and utilities can be modified by using the construction
`type_name_category`, where type is either "qaly" or "cost", name is the
name (e.g., "default") and category is the category used (e.g.,
"instant"), so one could pass `cost_default_instant` and modify the
cost. This will overwrite the value defined in the corresponding
cost/utility section.

The function is different from modify_item in that this function
evaluates sequentially the arguments within the list passed. This
implies a slower performance relative to modify_item, but it can be more
cleaner and convenient in certain instances.

This function is intended to be used only within the `add_reactevt`
function in its `input` parameter and should not be run elsewhere or it
will return an error.

## Examples

``` r
add_reactevt(name_evt = "idfs",input = {
  modify_item_seq(list(cost.idfs = 500, cost.tx = cost.idfs + 4000))
  })
#> $idfs
#> $idfs$react
#> {
#>     modify_item_seq(list(cost.idfs = 500, cost.tx = cost.idfs + 
#>         4000))
#> }
#> 
#> 
```
