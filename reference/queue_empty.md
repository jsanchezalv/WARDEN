# Check if the event queue is empty

Check if the event queue is empty

## Usage

``` r
queue_empty(ptr, exclude_inf = FALSE)
```

## Arguments

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- exclude_inf:

  Logical, whether to exclude events with Inf time. Default is FALSE.

## Value

Logical, TRUE if the queue is empty, FALSE otherwise.
