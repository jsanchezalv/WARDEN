# Add cost to list --------------------------------------------------------
#' Defining costs for events and intervention
#'
#' @param .data Existing cost data
#' @param cost Value or expression to calculate the cost estimate
#' @param evt Vector of events for which this cost is applicable
#' @param trt Vector of interventions for which this cost is applicable
#' @param cycle_l Cycle length; only needed if costs are calculated per cycle
#' @param cycle_starttime Cycle when costs start being accrued; only needed if costs are calculated per cycle
#'
#' @return A list of costs
#' @export
#'
#' @details
#' Costs can be defined by writing expressions and objects in the cost argument whose execution will be delayed until the model runs.
#'
#' This function accepts the use of pipes (%>%) to define multiple costs.
#'
#' @examples
#' add_cost(evt = c("start","idfs","ttot"),trt = "int",cost = cost.int*fl.int + cost.idfs)
#'
add_cost <- function(.data=NULL,cost,evt,trt,category="default",cycle_l=NULL,cycle_starttime=0){
  
  if (category=="default") {
    warning("No cost category name assigned, used 'default'")
  }
  category <- paste0("cost_",category)
  
  data_list <- .data
  for (trt_item in trt) {
    for (event_item in evt) {

      cost_item <- list(list(cost=substitute(cost),
                             cycle_l = substitute(cycle_l),
                             cycle_starttime = substitute(cycle_starttime),
                             category=category))
      names(cost_item) <- paste(event_item,trt_item,category,sep="_")

      if (is.null(data_list)) {
        data_list <- cost_item
      } else{
        data_list <- append(data_list,cost_item)
      }
    }
  }
  return(data_list)
}


# Add utility to list --------------------------------------------------------

#' Defining utilities for events and interventions
#'
#' @param .data Existing utility data
#' @param util Value or expression to calculate the utility estimate
#' @param evt Events for which this utility is applicable
#' @param trt Interventions for which this utility is applicable
#' @param cycle_l Cycle length; only needed if utilities are calculated per cycle
#' @param cycle_starttime Cycle when utilities start being accrued; only needed if utilities are calculated per cycle
#'
#' @return A list of utilities
#' @export
#'
#' @details
#' Utilities can be defined by writing expressions and objects in the cost argument whose execution will be delayed until the model runs.
#'
#' This function accepts the use of pipes (%>%) to define multiple utilities.
#'
#' @examples
#' add_util(evt = c("start","idfs","ttot"),
#' trt = c("int", "noint"),
#' util = util.idfs.ontx * fl.idfs.ontx + util.idfs.offtx * (1-fl.idfs.ontx))

add_util <- function(.data=NULL,util,evt,trt,category="default",cycle_l=NULL,cycle_starttime=0){

  if (category=="default") {
    warning("No utility category name assigned, used 'default'")
  }
  category <- paste0("qaly_",category)
  data_list <- .data
  for (trt_item in trt) {

    for (event_item in evt) {

      util_item <- list(list(util=substitute(util),
                             cycle_l = substitute(cycle_l),
                             cycle_starttime = substitute(cycle_starttime),
                             category = category))
      names(util_item) <- paste(event_item,trt_item,category,sep="_")

      if (is.null(data_list)) {
        data_list <- util_item
      } else{
        data_list <- append(data_list,util_item)
      }
    }
  }
  return(data_list)
}


# Add item/parameter to list --------------------------------------------------------

#' Defining parameters that may be used in model calculations
#'
#' @param .data Existing data
#' @param ... Items to define for the simulation
#'
#' @return A list of items
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `add_item` with a list of one element, it's better to group them into a single `add_item` with a list of two elements.
#'
#' @examples
#' add_item(fl.idfs = 0)
#' add_item(util_idfs = if(psa_bool){rnorm(1,0.8,0.2)} else{0.8}, util.mbc = 0.6, cost_idfs = 2500)
#'

add_item <- function(.data=NULL,...){

  data_list <- .data
  
  list_item <- as.list(substitute(...()))
  
  if (is.null(data_list)) {
    data_list <- list_item
  } else{
    data_list <- append(data_list,list_item)
  }
  
  
  return(data_list)
}



# Add event to list of events ---------------------------------------------

#' Generate new events to be added to existing vector of events
#'
#' @param evt Event name and event time
#' @param env_ch Environment in which to save list (should not be defined by user)
#'
#' @importFrom stats setNames
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `new_event` with a list of one element, it's better to group them into a single `new_event` with a list of two elements.
#'
#' @examples
#' \dontrun{
#' new_event(list("ae"=5))
#' new_event(list("ae"=5,"nat.death" = 100))
#' }

new_event <- function(evt, env_ch = NULL){
  new_evt_name <- names(evt)
  new_evt <- setNames(unlist(evt),new_evt_name)
  
  input_list_trt <- parent.frame()$input_list_trt
  
  evtlist_temp <- list(cur_evtlist = c(input_list_trt$cur_evtlist,
                                  new_evt))
  
  evtlist_temp$cur_evtlist <- sort(evtlist_temp$cur_evtlist) 
  
  input_list_trt[["cur_evtlist"]] <- evtlist_temp$cur_evtlist
  
  list2env(input_list_trt["cur_evtlist"],envir = parent.frame())
  assign("input_list_trt",input_list_trt, envir = parent.frame())

}

# Modify event in list of events ---------------------------------------------

#' Modify the time of existing events
#'
#' @param evt A list of events and their times
#' @param env_ch Environment in which to save list (should not be defined by user)
#'
#' @importFrom utils modifyList
#' @importFrom stats setNames
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_event` with a list of one element, it's better to group them into a single `modify_event` with a list of two elements.
#'
#'
#' @examples
#' \dontrun{
#' modify_event(list("os"=40, "ttot"=curtime+0.0001))
#' }

modify_event <- function(evt, env_ch = NULL){
  input_list_trt <- parent.frame()$input_list_trt
  
  names_obj_temp <- names(input_list_trt$cur_evtlist)
  names_evt <- names(evt)
  names_found <- names_evt[names_evt %in% names_obj_temp]
  if (length(names_found)==0) {
    warning("Some or all event/s in modify_evt within ", paste(names_evt,collapse=", "), " not found. Use new_evt to add new events.")
  }
  matched <- which(names_obj_temp %in% names_found)
  
  input_list_trt[["cur_evtlist"]][matched] <- unlist(evt[names_obj_temp[names_obj_temp %in% names_found]])
  
  input_list_trt[["cur_evtlist"]] <- sort(input_list_trt[["cur_evtlist"]])
  
  list2env(input_list_trt["cur_evtlist"],envir = parent.frame())
  assign("input_list_trt",input_list_trt, envir = parent.frame())
  

}





# Modify item in input list -------------------------------------------------------------------------------------------------------------------------------

#' Modify the value of existing items
#'
#' @param list_item A list of items and their values or expressions
#' @param env_ch Environment in which to save list (should not be defined by user)
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_item` with a list of one element, it's better to group them into a single `modify_item` with a list of two elements.
#'
#' @examples
#' \dontrun{
#' modify_item(list(cost.idfs = 500, cost.tx = cost.tx + 4000))
#' }

modify_item <- function(list_item, env_ch = NULL){
  input_list_trt <- parent.frame()$input_list_trt
  
  input_list_trt[names(list_item)] <- lapply(list_item, unname)
  
  list2env(list_item,envir = parent.frame())
  assign("input_list_trt",input_list_trt, envir = parent.frame())
}



# Add_reactevt -------------------------------------------------------------------------------------------------------------------------------------------

#' Define the modifications to other events, costs, utilities, or other items affected by the occurrence of the event
#'
#' @param .data Existing data for event reactions
#' @param name_evt Name of the event for which reactions are defined.
#' @param input Expressions that define what happens at the event, using functions as defined in the Details section
#'
#' @export
#'
#' @details
#' There are a series of objects that can be used in this context to help define the event reactions.
#'
#' The following functions may be used to define event reactions within this `add_reactevt()` function:
#' `modify_item()` | Adds & Modifies items/flags/variables for future events
#' `new_event()` | Adds events to the vector of events for that patient
#' `modify_event()` | Modifies existing events by changing their time
#'
#' Apart from the items defined with add_item(), we can also use standard variables that are always defined within the simulation:
#' `curtime` | Current event time (numeric)
#' `prevtime` | Time of the previous event (numeric)
#' `cur_evtlist` | Named vector of events that is yet to happen for that patient (named numeric vector)
#' `evt` | Current event being processed (character)
#' `i` | Patient being iterated (character)
#' `simulation` | Simulation being iterated (numeric)
#'
#' The model will run until `curtime` is set to `Inf`, so the event that terminates the model should modify `curtime` and set it to `Inf`.
#'
#' @examples
#' add_reactevt(name_evt = "start",input = {})
#' add_reactevt(name_evt = "idfs",input = {modify_item(list("fl.idfs"= 0))})

add_reactevt <- function(.data=NULL,name_evt,input){


  data_list <- .data

  evt_r <- list(list(react=substitute(input)))
  
  names(evt_r) <- paste(name_evt)

  if (is.null(data_list)) {
    data_list <- evt_r
  } else{
    data_list <- append(data_list,evt_r)
  }


  return(data_list)
}



# Add drawing and initial event list -------------------------------------------------------------------------------------------------------------------------------------------

#' Define events and the initial event time
#'
#' @param .data Existing data for initial event times
#' @param trt The intervention for which the events and initial event times are defined
#' @param evts A vector of the names of the events
#' @param other_inp A vector of other input variables that should be saved during the simulation
#' @param input The definition of initial event times for the events listed in the evts argument
#'
#' @return A list of initial events and event times
#' @export
#'
#' @details
#' Events need to be separately defined for each intervention.
#'
#' For each event that is defined in this list, the user needs to add a reaction to the event using the `add_reactevt()` function which will determine what calculations will happen at an event.
#'
#' @examples
#' add_tte(trt="int",evts = c("start","ttot","idfs","os"),
#' input={
#' start <- 0
#' idfs <- draw_tte(1,'lnorm',coef1=2, coef2=0.5)
#' ttot <- min(draw_tte(1,'lnorm',coef1=1, coef2=4),idfs)
#' os <- draw_tte(1,'lnorm',coef1=0.8, coef2=0.2)
#' })
#'
add_tte <- function(.data=NULL,trt, evts, other_inp = NULL,input){
  data_list <- .data

  for (trt in trt) {

    if (!is.character(other_inp) & !is.null(other_inp)) {
      stop("other_inp argument is required to be a character vector or be set to NULL")
    }

    if (!is.character(evts) | length(evts)<2) {
      stop("evts argument in add_tte for the intervention ", trt, " is required to be a character vector with length >1")
    }


    evt_l <- list(list(expr=substitute(input),
                       evts = evts,
                       other_inp = other_inp
    ))
    names(evt_l) <- paste(trt)

    if (is.null(data_list)) {
      data_list <- evt_l
    } else{
      data_list <- append(data_list,evt_l)
    }

  }

  if (length(names(data_list))!=length(unique(names(data_list)))) {
    stop("More than one initial set of time to events has been defined for an intervention.")
  }

  return(data_list)
}



# Continuous and instantaneous discounting ----------------------------------------------------------------------------------------------------------------

#' Calculate discounted costs and qalys between events
#'
#' @param lcldr The discount rate 
#' @param lclprvtime The time of the previous event in the simulation
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on continuous time discounting
#'
#' @examples AddOngoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#'
#' @export
#' @details

AddOngoing <- function(lcldr=0.035, lclprvtime, lclcurtime, lclval){
  
  Instantdr <- log(1+lcldr)
  
  # calculate additional qalys
  if (lclprvtime==lclcurtime | is.null(lclval)) {
    add <- 0
  } else if(lcldr==0) {
    add <- lclval*(lclcurtime - lclprvtime)
  } else{
    add <- ((lclval)/(0 - Instantdr)) * (exp(lclcurtime * ( 0 - Instantdr)) - exp(lclprvtime * (0 - Instantdr)))
    
  }
  
  return(add)
}


#' Calculate instantaneous discounted costs or qalys
#'
#' @param lcldr The discount rate
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on discrete time discounting
#'
#' @examples AddInstant(lcldr=0.035, lclcurtime=3, lclval=2500)
#'
#' @export
#' 
AddInstant <- function(lcldr=0.035, lclcurtime, lclval){
  
  if (is.null(lclval)) {
    addinst <- 0
  } else{
    addinst <- lclval * ((1+lcldr)^(-lclcurtime))    
  }# Note use of DISCRETE TIME discounting for instantaneous costs and benefits
  
  return(addinst)
}


# Cycle discounting -------------------------------------------------------

#' Cycle discounting
#'
#' @param lcldr The discount rate
#' @param lclprvtime The time of the previous event in the simulation
#' @param cyclelength The cycle length
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The  value to be discounted
#' @param starttime The start time for accrual of cycle costs (if not 0)
#'
#' @return Double based on cycle discounting
#'
#' @examples AddCycle(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#'
#' @export

AddCycle <- function(lcldr=0.035, lclprvtime=0, cyclelength,lclcurtime, lclval,starttime=0){
  
  addcycle <- 0
  
  if (is.null(lclval)) {
    addcycle <- 0
  } else{
    
    #Note this makes the cycle utilities work weird, so do not use cycle utilities for now!
    for (i in 1:length(lclval)) {
      lclval_i <- lclval[i]
      starttime_i <- starttime[i]
      cyclelength_i <- cyclelength[i]
      
      
      if (lclval_i==0 ) {} else{
        
        cycle.time.total <- if(starttime_i>= lclcurtime){0}else{seq(from=starttime_i, to = lclcurtime , by= cyclelength_i)} #all cycles that happened until current time of event
        
        # cycle.time.total <- seq(from=starttime, to = lclcurtime , by= cyclelength) #all cycles that happened until current time of event
        
        #If the cost starts at the selected starttime or at time 0, then include that time, otherwise exclude it
        if (lclprvtime==0) {
          cycle.time <- c(0,cycle.time.total[cycle.time.total >= lclprvtime])  #times at which the cycles take place during this event, put this condition to count also time 0
          n_cycles <- length(cycle.time)
          s <- (1+lcldr)^cyclelength_i -1
          addcycle <- sum(addcycle,lclval_i * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1) )
          
        } else{
          if (starttime_i ==lclprvtime) {
            cycle.time <- cycle.time.total[cycle.time.total >= lclprvtime]  #times at which the cycles take place during this event, put this condition to count also time of the previous event
          } else{
            cycle.time <- cycle.time.total[cycle.time.total > lclprvtime]  #times at which the cycles take place during this event
          }
          n_cycles_remaining <- length(cycle.time)
          d <- lclprvtime/cyclelength_i
          s <- (1+lcldr)^cyclelength_i -1
          addcycle <- sum(addcycle, lclval_i * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d)) )
        }
        
        #If starting from 0, can be changed substituting interest rate such that s = (1+r)^cyclelength - 1, and using the formula that lclvalq * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1)
        #If starting from time t, then compute transformed time as d = t/cyclelength and use lclvalq * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d-1)), where
        #n_cycles_remaining is the n_cycles - d (so the remaining cycles to be considered), e.g. if 13 cycles (From t=0), and delay 6 periods, then n_cycles_remaining = 7 and d=6
        
        # combine additional costs and additional quantity in a list
        
      }
      
    }
  }
  return(addcycle)
}