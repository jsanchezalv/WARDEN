# Calculate instantaneous discounted costs or qalys for vectors

Calculate instantaneous discounted costs or qalys for vectors

## Usage

``` r
disc_instant_v(lcldr, lclcurtime, lclval)
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
disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500)
#> [1] 2254.857

```
