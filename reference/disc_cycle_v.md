# Cycle discounting for vectors

Cycle discounting for vectors

## Usage

``` r
disc_cycle_v(
  lcldr,
  lclprvtime,
  cyclelength,
  lclcurtime,
  lclval,
  starttime,
  max_cycles = NULL
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

- max_cycles:

  The maximum number of cycles

## Value

Double vector based on cycle discounting

## Details

This function per cycle discounting, i.e., considers that the cost/qaly
is accrued per cycles, and performs it automatically without needing to
create new events. It can accommodate changes in cycle
length/value/starttime (e.g., in the case of induction and maintenance
doses) within the same item.

## Examples

``` r
disc_cycle_v(lcldr=0.03, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#> [1] 11666.54
disc_cycle_v(
 lcldr=0.000001,
 lclprvtime=0,
 cyclelength=1/12,
 lclcurtime=2,
 lclval=500,
 starttime=0,
 max_cycles = 4)
#> [1] 2000

#Here we have a change in cycle length, max number of cylces and starttime at time 2
 #(e.g., induction to maintenance)
#In the model, one would do this by redifining cycle_l, max_cycles and starttime
 #of the corresponding item at a given event time. 
disc_cycle_v(lcldr=0,
 lclprvtime=c(0,1,2,2.5),
 cyclelength=c(1/12, 1/12,1/2,1/2),
 lclcurtime=c(1,2,2.5,4), lclval=c(500,500,500,500),
 starttime=c(0,0,2,2), max_cycles = c(24,24,2,2)
  )
#> [1] 6000 6000  500  500
```
