# Get a specific event time

Get a specific event time

## Usage

``` r
get_event(event_name, ptr, patient_id)
```

## Arguments

- event_name:

  Character string, the name of the event.

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- patient_id:

  The patient ID. Defaults to `i`.

## Value

Numeric, time of event for patient
