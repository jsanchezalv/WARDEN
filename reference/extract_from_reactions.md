# Extract all items and events and their interactions from the event reactions list

Extract all items and events and their interactions from the event
reactions list

## Usage

``` r
extract_from_reactions(reactions)
```

## Arguments

- reactions:

  list generated through `add_reactevt`

## Value

A data.frame with the relevant item/event, the event where it's
assigned, and whether it's contained within a conditional statement

## Examples

``` r
evt_react_list2 <-
  add_reactevt(name_evt = "sick",
               input = {modify_item(list(a=1+5/3))
                 assign("W", 5 + 3 / 6 )
                 x[5] <- 18
                 for(i in 1:5){
                   assign(paste0("x_",i),5+3)
                 }
                 if(j == TRUE){
                   y[["w"]] <- 612-31+3
                 }#'                
                 q_default <- 0
                 c_default <- 0
                 curtime   <- Inf
                 d <- c <- k <- 67
               })
   
 extract_from_reactions(evt_react_list2)
#>      event            name   type conditional_flag   definition
#>     <char>          <char> <char>           <lgcl>       <char>
#>  1:   sick               a   item            FALSE      1 + 5/3
#>  2:   sick               W   item            FALSE      5 + 3/6
#>  3:   sick            x[5]   item            FALSE           18
#>  4:   sick paste0('x_', i)   item            FALSE        5 + 3
#>  5:   sick        y[['w']]   item             TRUE 612 - 31 + 3
#>  6:   sick       q_default   item            FALSE            0
#>  7:   sick       c_default   item            FALSE            0
#>  8:   sick         curtime   item            FALSE          Inf
#>  9:   sick               d   item            FALSE c <- k <- 67
#> 10:   sick               c   item            FALSE      k <- 67
#> 11:   sick               k   item            FALSE           67

```
