# Release a discrete resource

Frees the resource for the current patient (`i`) and, if the patient was
actually using the resource and a `resume_event` is supplied, schedules
that event for the next patient in the queue.

## Usage

``` r
release(resource, resume_event = NULL, amount = 1L)
```

## Arguments

- resource:

  A `resource_discrete` object.

- resume_event:

  Character string. Name of the event to schedule for the next queued
  patient, or `NULL` (default) to skip re-triggering.

- amount:

  Integer. Units to free (default `1L`).

## Value

Invisibly, `TRUE` if the patient was using the resource, `FALSE`
otherwise.
