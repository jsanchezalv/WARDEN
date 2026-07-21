# Scale remaining time to an event by a hazard ratio

Computes `curtime + (event_time - curtime) / hr`. Use inside event
reactions to shrink or expand the remaining time to an event.

## Usage

``` r
scale_remaining_time(event_name, hr)
```

## Arguments

- event_name:

  Character. Name of the event whose time to scale.

- hr:

  Numeric greater than 0. Hazard ratio to apply. Values below 1 extend
  time; values above 1 shorten it.

## Value

Numeric: the new absolute event time.
