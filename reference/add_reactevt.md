# Define the modifications to other events, costs, utilities, or other items affected by the occurrence of the event

Define the modifications to other events, costs, utilities, or other
items affected by the occurrence of the event

## Usage

``` r
add_reactevt(.data = NULL, name_evt, input)
```

## Arguments

- .data:

  Existing data for event reactions

- name_evt:

  Name of the event for which reactions are defined.

- input:

  Expressions that define what happens at the event, using functions as
  defined in the Details section

## Value

A named list with the event name, and inside it the substituted
expression saved for later evaluation

## Details

There are a series of objects that can be used in this context to help
define the event reactions.

The following functions may be used to define event reactions within
this `add_reactevt()` function:
[`modify_item()`](https://jsanchezalv.github.io/WARDEN/reference/modify_item.md)
\| Adds & Modifies items/flags/variables for future events (does not
consider sequential)
[`modify_item_seq()`](https://jsanchezalv.github.io/WARDEN/reference/modify_item_seq.md)
\| Adds & Modifies items/flags/variables for future events in a
sequential manner
[`new_event()`](https://jsanchezalv.github.io/WARDEN/reference/new_event.md)
\| Adds events to the vector of events for that patient
[`modify_event()`](https://jsanchezalv.github.io/WARDEN/reference/modify_event.md)
\| Modifies existing events by changing their time

Apart from the items defined with add_item(), we can also use standard
variables that are always defined within the simulation: `curtime` \|
Current event time (numeric) `prevtime` \| Time of the previous event
(numeric) `cur_evtlist` \| Named vector of events that is yet to happen
for that patient (named numeric vector) `evt` \| Current event being
processed (character) `i` \| Patient being iterated (character)
`simulation` \| Simulation being iterated (numeric)

The model will run until `curtime` is set to `Inf`, so the event that
terminates the model should modify `curtime` and set it to `Inf`.

The user can use `extract_from_reactions` function on the output to
obtain a data.frame with all the relationships defined in the reactions
in the model.

## Examples

``` r
add_reactevt(name_evt = "start",input = {})
#> $start
#> $start$react
#> {
#> }
#> 
#> 
add_reactevt(name_evt = "idfs",input = {modify_item(list("fl.idfs"= 0))})
#> $idfs
#> $idfs$react
#> {
#>     modify_item(list(fl.idfs = 0))
#> }
#> 
#> 
```
