# Calculate discounted costs and qalys between events

Calculate discounted costs and qalys between events

## Usage

``` r
disc_ongoing(lcldr = 0.035, lclprvtime, lclcurtime, lclval)
```

## Arguments

- lcldr:

  The discount rate

- lclprvtime:

  The time of the previous event in the simulation

- lclcurtime:

  The time of the current event in the simulation

- lclval:

  The value to be discounted

## Value

Double based on continuous time discounting

## Examples

``` r
disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#> Warning: `disc_ongoing()` was deprecated in WARDEN 2.0.0.
#> ℹ Please use `disc_ongoing_v()` instead.
#> [1] 5886.65

```
