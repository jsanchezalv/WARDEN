# Release multiple discrete resources

Frees the current patient (`i`) from each resource via a single C++
call. Removes the patient from both the using list and the queue (all
entries). Schedules resume events only for resources where the patient
was actually using (i.e., where capacity was freed).

## Usage

``` r
release_all(resources, resume_event = NULL, amounts = NULL)
```

## Arguments

- resources:

  A list of `resource_discrete` objects.

- resume_event:

  `NULL` (no re-triggering), a single string (same event for all
  resources' queues), or a character vector of length
  `length(resources)` (per-resource events, matched positionally).

- amounts:

  Integer vector of units per resource (default `1L` for each).

## Value

Invisibly `NULL`.
