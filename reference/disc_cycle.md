# Cycle discounting

Cycle discounting

## Usage

``` r
disc_cycle(
  lcldr = 0.035,
  lclprvtime = 0,
  cyclelength,
  lclcurtime,
  lclval,
  starttime = 0
)
```

## Arguments

- lcldr:

  The discount rate

- lclprvtime:

  The time of the previous event in the simulation

- cyclelength:

  The cycle length

- lclcurtime:

  The time of the current event in the simulation

- lclval:

  The value to be discounted

- starttime:

  The start time for accrual of cycle costs (if not 0)

## Value

Double based on cycle discounting

## Details

Note this function counts both extremes of the interval, so the example
below would consider 25 cycles, while disc_cycle_v leave the right
interval open

## Examples

``` r
disc_cycle(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#> Warning: `disc_cycle()` was deprecated in WARDEN 2.0.0.
#> ℹ Please use `disc_cycle_v()` instead.
#> [1] 12079.88

```
