# Pop and Return the Next Event

Removes the next event from the queue and returns its details. Not
needed by user.

## Usage

``` r
pop_and_return_event(ptr)
```

## Arguments

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

## Value

A named list with `patient_id`, `event_name`, and `time`.
