# Decrement a shared counter

Subtracts `delta` from the value of a `shared_input` object and returns
the new value.

## Usage

``` r
shared_decr(obj, delta = 1)
```

## Arguments

- obj:

  A `shared_input` object (constrained / shared mode).

- delta:

  Numeric. Amount to subtract (default `1`).

## Value

The new value after decrementing.
