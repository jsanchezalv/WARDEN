# Remove Events for a Patient

Removes one or more events from the queue for the given patient.

## Usage

``` r
remove_event(events, ptr, patient_id)
```

## Arguments

- events:

  A character vector of event names to remove. It can also handle lists
  instead of named vectors (at a small computational cost).

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- patient_id:

  The patient ID. Defaults to `i`.

## Value

NULL (invisible). Modifies the queue in-place.
