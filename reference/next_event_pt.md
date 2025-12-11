# Get the Next Event(s) in the Queue for a specific patient

Retrieves the next `n` events (without removing them).

## Usage

``` r
next_event_pt(n = 1, ptr, patient_id)
```

## Arguments

- n:

  Number of events to retrieve. Default is 1.

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- patient_id:

  The patient ID. Defaults to `i`.

## Value

A list of events, each with `patient_id`, `event_name`, and `time`.
