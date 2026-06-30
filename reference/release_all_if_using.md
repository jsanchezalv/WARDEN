# Release multiple discrete resources (using only)

Frees the current patient (`i`) from each resource only if they are
currently using it. Applies an all_or_none policy: only acts if the
patient is using ALL listed resources. Does not remove the patient from
any queue. Use this when a patient transitions state but should keep
their queue position.

## Usage

``` r
release_all_if_using(resources, resume_event = NULL, amounts = NULL)
```

## Arguments

- resources:

  A list of `resource_discrete` objects.

- resume_event:

  `NULL` (no re-triggering), a single string (same event for all
  resources' queues), or a character vector of length
  `length(resources)` (per-resource events, matched positionally).

- amounts:

  Integer vector of units per resource, or `NULL` (default). `NULL`
  releases all entries and purges queue entries for the patient.
  Specified amounts perform exact indivisible releases without purging.

## Value

Invisibly `NULL`.
