# Compute effective hazard ratio with treatment waning

Returns the effective hazard ratio at a given time on treatment,
accounting for waning of treatment effect.

## Usage

``` r
waning_hr(base_hr, time_on_tx, start, duration, type = "linear")
```

## Arguments

- base_hr:

  Numeric. Initial treatment HR (\< 1 means benefit).

- time_on_tx:

  Numeric. Time the patient has been on treatment.

- start:

  Numeric. Time at which waning begins (HR = base_hr before this).

- duration:

  Numeric. Duration of the waning period (HR reaches 1 at start +
  duration).

- type:

  Character. Waning function: `"linear"` (default), `"exponential"`, or
  `"instant"`.

## Value

Numeric scalar: effective HR at `time_on_tx`.
