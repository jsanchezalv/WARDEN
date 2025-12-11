# Add Events to the Queue for a Patient

Adds one or more events for a given patient to the queue.

## Usage

``` r
new_event(events, ptr, patient_id)
```

## Arguments

- events:

  A named numeric vector. Names are event types, values are event times.
  It can also handle lists instead of named vectors (at a small
  computational cost).

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- patient_id:

  The patient ID. Defaults to `i`.

## Value

NULL (invisible). Modifies the queue in-place.

## Details

The functions to add/modify events/inputs use named vectors or lists.
Whenever several inputs/events are added or modified, it's recommended
to group them within one function, as it reduces the computation cost.
So rather than use two `new_event` with a list of one element, it's
better to group them into a single `new_event` with a list of two
elements.

While multiple events can be added, they must be named differently. If
the same event is added multiple times at once, only the last occurrence
will be kept (only one event per event type in the queue of events yet
to occur). If an event occurs, then a new one with the same name can be
set.

This function is intended to be used only within the `add_reactevt`
function in its `input` parameter and should not be run elsewhere or it
will return an error.

## Examples

``` r
add_reactevt(name_evt = "idfs",input = {new_event(c("ae"=5))})
#> $idfs
#> $idfs$react
#> {
#>     new_event(c(ae = 5))
#> }
#> 
#> 
```
