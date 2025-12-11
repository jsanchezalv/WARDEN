# Check if a Patient Has a Specific Event

Check if a Patient Has a Specific Event

## Usage

``` r
has_event(event_name, ptr, patient_id, exclude_inf = FALSE)
```

## Arguments

- event_name:

  Character string, the name of the event.

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- patient_id:

  The patient ID. Defaults to `i`.

- exclude_inf:

  Logical, whether to exclude events with Inf time. Default is FALSE.

## Value

Logical, TRUE if the event exists for the patient (optionally excluding
Inf), FALSE otherwise.
