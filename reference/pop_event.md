# Remove the Next Event from the Queue

Removes the next scheduled event from the queue. Not needed by user.

## Usage

``` r
pop_event(ptr)
```

## Arguments

- ptr:

  The event queue pointer. Defaults to `cur_evtlist`.

## Value

NULL (invisible). Modifies the queue in-place.
