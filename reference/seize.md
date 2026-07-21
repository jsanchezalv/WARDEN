# Seize a discrete resource

Convenience wrapper around `resource$attempt_block()`. Reads `i` and
`curtime` from the calling environment automatically.

## Usage

``` r
seize(resource, amount = 1L, priority = 1L, ...)
```

## Arguments

- resource:

  A `resource_discrete` object.

- amount:

  Integer. Number of resource units to seize (default `1L`).

- priority:

  Integer. Queue priority for the patient (default `1L`). Lower values
  are higher priority.

- ...:

  Not used. Passing
  [`seize_all()`](https://jsanchezalv.github.io/WARDEN/reference/seize_all.md)-specific
  arguments here (e.g., `accum_queue`, `force_unblock`, `policy`) will
  raise an error.

## Value

`TRUE` if acquired, `FALSE` if queued, `NA` if rejected (queue full,
only when `max_queue` is set on the resource).
