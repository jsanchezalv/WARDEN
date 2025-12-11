# Calculate instantaneous discounted costs or qalys

Calculate instantaneous discounted costs or qalys

## Usage

``` r
disc_instant(lcldr = 0.035, lclcurtime, lclval)
```

## Arguments

- lcldr:

  The discount rate

- lclcurtime:

  The time of the current event in the simulation

- lclval:

  The value to be discounted

## Value

Double based on discrete time discounting

## Examples

``` r
disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500)
#> Warning: `disc_instant()` was deprecated in WARDEN 2.0.0.
#> ℹ Please use `disc_instant_v()` instead.
#> [1] 2254.857

```
