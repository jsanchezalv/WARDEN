# Seize a discrete resource

Convenience wrapper around `resource$attempt_block()`. Reads `i` and
`curtime` from the calling environment automatically.

## Usage

``` r
seize(resource, amount = 1L)
```

## Arguments

- resource:

  A `resource_discrete` object.

- amount:

  Integer. Number of resource units to seize (default `1L`).

## Value

`TRUE` if acquired, `FALSE` if queued, `NA` if rejected (queue full,
only when `max_queue` is set on the resource).
