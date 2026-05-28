# Seize multiple discrete resources atomically

Attempts to acquire a list of resources in a single C++ call, avoiding
per-resource R-to-C++ round trips.

## Usage

``` r
seize_all(resources, policy = "all_or_none", amounts = NULL, priorities = NULL)
```

## Arguments

- resources:

  A list of `resource_discrete` objects.

- policy:

  `"all_or_none"` (default) — acquires all only if all are available;
  queues for the first bottleneck without holding any others.
  `"sequential"` — acquires in list order, holding partial acquisitions;
  user is responsible for avoiding deadlocks.

- amounts:

  Integer vector of units per resource (default `1L` for each).

- priorities:

  Integer vector of priorities per resource (default `1L`).

## Value

`TRUE` if all acquired, `FALSE` if queued for a bottleneck resource,
`NA` if rejected.
