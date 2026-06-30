# Release multiple discrete resources

Frees the current patient (`i`) from each resource. Applies an
all_or_none policy: only acts if the patient is using ALL listed
resources. Removes the patient from both the using list and the queue.
When `amounts` is `NULL` (default), releases all entries and purges
queue entries for the patient; when `amounts` is specified, performs an
exact indivisible release without purging the queue. Schedules resume
events only for resources where the patient was actually using (i.e.,
where capacity was freed).

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

  Integer vector of units per resource, or `NULL` (default). `NULL`
  releases all entries and purges queue entries for the patient.
  Specified amounts perform exact indivisible releases without purging.

## Value

Invisibly `NULL`.
