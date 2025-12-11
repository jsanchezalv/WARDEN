# Get the Next Event(s) in the Queue

Retrieves the next `n` events (without removing them).

## Usage

``` r
next_event(n = 1, ptr)
```

## Arguments

- n:

  Number of events to retrieve. Default is 1.

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

## Value

A list of events, each with `patient_id`, `event_name`, and `time`.
