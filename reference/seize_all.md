# Seize multiple discrete resources atomically

Attempts to acquire a list of resources in a single C++ call, avoiding
per-resource R-to-C++ round trips.

## Usage

``` r
seize_all(
  resources,
  policy = "all_or_none",
  amounts = NULL,
  priorities = NULL,
  force_unblock = FALSE,
  accum_queue = TRUE
)
```

## Arguments

- resources:

  A list of `resource_discrete` objects.

- policy:

  `"all_or_none"` (default) – acquires all only if all are available;
  queues for the first bottleneck without holding any others.
  `"sequential"` – acquires in list order, holding partial acquisitions;
  user is responsible for avoiding deadlocks.

- amounts:

  Integer vector of units per resource (default `1L` for each).

- priorities:

  Integer vector of priorities per resource (default `1L`).

- force_unblock:

  Logical (default `FALSE`). When `TRUE` and all resources have
  sufficient capacity but the patient is not first in queue on some
  resources (deadlock), forcibly moves the patient to the front of those
  queues and acquires all resources. Has no effect when capacity is
  genuinely insufficient.

- accum_queue:

  Logical (default `TRUE`). When `FALSE` and the patient already has a
  queue entry on a bottleneck resource from a previous failed attempt,
  skips adding a new entry (existing entry serves). Prevents phantom
  queue accumulation on retries.

## Value

`TRUE` if all acquired, `FALSE` if queued for a bottleneck resource,
`NA` if rejected.
