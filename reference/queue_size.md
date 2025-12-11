# Get the Size of the Event Queue

Get the Size of the Event Queue

## Usage

``` r
queue_size(ptr, exclude_inf = FALSE)
```

## Arguments

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

- exclude_inf:

  Logical, whether to exclude events with Inf time. Default is FALSE.

## Value

An integer indicating the number of events in the queue.
