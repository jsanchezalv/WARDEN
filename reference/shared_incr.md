# Increment a shared counter

Adds `delta` to the value of a `shared_input` object and returns the new
value.

## Usage

``` r
shared_incr(obj, delta = 1)
```

## Arguments

- obj:

  A `shared_input` object (constrained / shared mode).

- delta:

  Numeric. Amount to add (default `1`).

## Value

The new value after incrementing.
