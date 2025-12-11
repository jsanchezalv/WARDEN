# Calculate discounted costs and qalys between events for vectors

Calculate discounted costs and qalys between events for vectors

## Usage

``` r
disc_ongoing_v(lcldr, lclprvtime, lclcurtime, lclval)
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
disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#> [1] 5886.65

```
