# Release a discrete resource only if using

Releases the resource for the current patient (`i`) only if they are
currently using it. Does nothing if the patient is not using the
resource. Never removes queue entries.

## Usage

``` r
release_if_using(resource, resume_event = NULL, amount = NULL)
```

## Arguments

- resource:

  A `resource_discrete` object.

- resume_event:

  Character string. Name of the event to schedule for the next queued
  patient, or `NULL` (default) to skip re-triggering.

- amount:

  Integer or `NULL` (default). Units to release. `NULL` means 1. Must
  exactly match the amount used when seizing (indivisible units).

## Value

Invisibly `NULL`.
