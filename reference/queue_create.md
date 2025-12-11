# Create a New Event Queue

Initializes a new event queue with the specified priority order of event
names.

## Usage

``` r
queue_create(priority_order)
```

## Arguments

- priority_order:

  A character vector of event names sorted by decreasing importance.

## Value

An external pointer to the new event queue.
