# Define events and the initial event time

Define events and the initial event time

## Usage

``` r
add_tte(.data = NULL, arm, evts, other_inp = NULL, input)
```

## Arguments

- .data:

  Existing data for initial event times

- arm:

  The intervention for which the events and initial event times are
  defined

- evts:

  A vector of the names of the events

- other_inp:

  A vector of other input variables that should be saved during the
  simulation

- input:

  The definition of initial event times for the events listed in the
  evts argument

## Value

A list of initial events and event times

## Details

Events need to be separately defined for each intervention.

For each event that is defined in this list, the user needs to add a
reaction to the event using the
[`add_reactevt()`](https://jsanchezalv.github.io/WARDEN/reference/add_reactevt.md)
function which will determine what calculations will happen at an event.

## Examples

``` r
add_tte(arm="int",evts = c("start","ttot","idfs","os"),
input={
start <- 0
idfs <- draw_tte(1,'lnorm',coef1=2, coef2=0.5)
ttot <- min(draw_tte(1,'lnorm',coef1=1, coef2=4),idfs)
os <- draw_tte(1,'lnorm',coef1=0.8, coef2=0.2)
})
#> $int
#> $int$expr
#> {
#>     start <- 0
#>     idfs <- draw_tte(1, "lnorm", coef1 = 2, coef2 = 0.5)
#>     ttot <- min(draw_tte(1, "lnorm", coef1 = 1, coef2 = 4), idfs)
#>     os <- draw_tte(1, "lnorm", coef1 = 0.8, coef2 = 0.2)
#> }
#> 
#> $int$evts
#> [1] "start" "ttot"  "idfs"  "os"   
#> 
#> $int$other_inp
#> NULL
#> 
#> 
```
