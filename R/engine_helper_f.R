# Initial event list --------------------------------------------------------------------------------------------------------------------------------------

#' Execute the initial time to events and separate the events from other inputs that are stored
#'
#' @param trt_name A character string of the name of the intervention
#' @param input_list_trt A list of simulation inputs
#'
#' @return A named vector of initial event times, and a named vector of other inputs to be stored
#'
#' @examples
#' InitEventList(trt = "int",input_list_trt = input_list_trt)
#'
#' @keywords internal
#' @noRd


InitEventList <- function(trt_name,input_list_trt){
  position <- which(trt_name==names(input_list_trt$init_event_list))

  time_data <- local({
    evts_v <- input_list_trt$init_event_list[[position]][["evts"]]
    othert_v <- input_list_trt$init_event_list[[position]][["other_inp"]]
    list2env(mget(c(evts_v,othert_v),ifnotfound=Inf), envir=environment()) #initialize
    eval(input_list_trt$init_event_list[[position]][["expr"]]) #run script
    evttime <- lapply(mget(evts_v,ifnotfound=Inf),unname) #get event times and make sure they are unnamed
    othertime <- if(!is.null(othert_v)){mget(othert_v,ifnotfound=Inf)} else{NULL}  #get other inputs times
    out <- list(evttime=evttime, othertime=othertime)
  },input_list_trt)

  #Event data
  cur_evtlist <- unlist(time_data$evttime)
  cur_evtlist <- sort(cur_evtlist)
  return(list(cur_evtlist = cur_evtlist, time_data = unlist(time_data$othertime)))
}



# Get next event ------------------------------------------------------------------------------------------------------------------------------------------

#' Identify which event to process next from a list of events
#'
#' @param evt_list A list of possible events with event times
#'
#' @return Two lists: one containing the name and time of the next event, the other with the remaining events to be processed
#'
#' @examples
#' GetNxtEvt(evt_list = input_list_trt$cur_evtlist)
#'
#' @keywords internal
#' @noRd

GetNxtEvt <- function(evt_list){                  # This function identifies which event is to be processed next for each patient, depending on intervention

  if (length(evt_list)>0) {
    cur_evtlist <- list(out = list(evt = names(evt_list[1]), evttime = evt_list[[1]]), evt_list = evt_list[-1])
  } else {
    cur_evtlist <- NULL
  }


  return(cur_evtlist)
}



# Reaction to Event ---------------------------------------------------------------------------------------------------------------------------------------

#' Evaluates the reactions of the event identified by GetNextEvt
#'
#' @param thisevt A two element list containing the first list from GetNextEvt: evt and evttime
#' @param trt A character string of the name of the intervention currently being processed
#' @param input_list_trt A list of simulation inputs
#'
#' @return The updated input list with after the reaction to the event is evaluated
#'
#' @examples
#' ReactEvt(thisevt="evt1",trt="int",input_list_trt=input_list_trt)
#'
#' @keywords internal
#' @noRd

ReactEvt <- function(thisevt,trt,input_list_trt=NULL){      # This function processes the next event (as identified in the GetNextEvt function)
  # Initial set-up --------------------------
  evt <- thisevt$evt                  # Identify event type
  prevtime <- input_list_trt$curtime                 # Identify time of previous event
  curtime <- thisevt$evttime         # Identify time of next event
  
  input_list_trt[["curtime"]] <- curtime
  input_list_trt[["evt"]] <- evt
  input_list_trt[["prevtime"]] <- prevtime
  input_list_trt[["trt"]] <- trt
  
  # Create costs and utilities for event --------------------------------------------------
  evt_trt <- paste(evt,trt,sep = "_")
  
  #For each cost/utility category, evaluate the relevant equation and get the raw inputs for ongoing/instant/cycle
  # Costs -------------------------------------------------------------------
  for (cost_cat in input_list_trt$uc_lists$cost_categories_ongoing) {
    input_list_trt[paste0(cost_cat,"_","ongoing")] <- get_input(input_list_trt$uc_lists$cost_ongoing_list,
                                                                ifnull=0,
                                                                type="cost",
                                                                evt_trt_i=paste(evt_trt,cost_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
  }
  for (cost_cat in input_list_trt$uc_lists$cost_categories_instant) {
    input_list_trt[paste0(cost_cat,"_","instant")] <- get_input(input_list_trt$uc_lists$cost_instant_list,
                                                                ifnull=0,
                                                                type="cost",
                                                                evt_trt_i=paste(evt_trt,cost_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
  }
  for (cost_cat in input_list_trt$uc_lists$cost_categories_cycle) {
    input_list_trt[paste0(cost_cat,"_","cycle")] <- get_input(input_list_trt$uc_lists$cost_cycle_list,
                                                              ifnull=0,
                                                              type="cost",evt_trt_i=paste(evt_trt,cost_cat,sep="_"),
                                                              input_list_trt_i=input_list_trt)
    input_list_trt[paste0(cost_cat,"_","cycle_l")] <- get_input(input_list_trt$uc_lists$cost_cycle_list,
                                                                ifnull=1,
                                                                type="cycle_l",
                                                                evt_trt_i=paste(evt_trt,cost_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
    input_list_trt[paste0(cost_cat,"_","cycle_starttime")] <- get_input(input_list_trt$uc_lists$cost_cycle_list,
                                                                        ifnull=0,
                                                                        type="cycle_starttime",
                                                                        evt_trt_i=paste(evt_trt,cost_cat,sep="_"),
                                                                        st_trt_i=input_list_trt)
  }
  
  # Utilities -------------------------------------------------------------------
  for (util_cat in input_list_trt$uc_lists$util_categories_ongoing) {
    input_list_trt[paste0(util_cat,"_","ongoing")] <- get_input(input_list_trt$uc_lists$util_ongoing_list,
                                                                ifnull=0,
                                                                type="util",
                                                                evt_trt_i=paste(evt_trt,util_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
  }
  for (util_cat in input_list_trt$uc_lists$util_categories_instant) {
    input_list_trt[paste0(util_cat,"_","instant")] <- get_input(input_list_trt$uc_lists$util_instant_list,
                                                                ifnull=0,
                                                                type="util",evt_trt_i=paste(evt_trt,util_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
  }
  for (util_cat in input_list_trt$uc_lists$util_categories_cycle) {
    input_list_trt[paste0(util_cat,"_","cycle")] <- get_input(input_list_trt$uc_lists$util_cycle_list,ifnull=0,
                                                              type="util",
                                                              evt_trt_i=paste(evt_trt,util_cat,sep="_"),
                                                              input_list_trt_i=input_list_trt)
    input_list_trt[paste0(util_cat,"_","cycle_l")] <- get_input(input_list_trt$uc_lists$util_cycle_list,
                                                                ifnull=1,
                                                                type="cycle_l",
                                                                evt_trt_i=paste(evt_trt,util_cat,sep="_"),
                                                                input_list_trt_i=input_list_trt)
    input_list_trt[paste0(util_cat,"_","cycle_starttime")] <- get_input(input_list_trt$uc_lists$util_cycle_list,
                                                                        ifnull=0,
                                                                        type="cycle_starttime",
                                                                        evt_trt_i=paste(evt_trt,util_cat,sep="_"),
                                                                        sinput_list_trt_i=input_list_trt)
  }

  #Evaluate the reaction to the event
  input_list_trt <- eval_reactevt(input_list_trt$evt_react_list, evt,input_list_trt)
  

  return(input_list_trt)

}


# Evaluate event ------------------------------------------------------------------------------------------------------------------------------------------


#' Calculates the expression which has been defined in the reaction of the event and computes discounting
#'
#' @param x The evt_react_list from the input_list_trt object. It contains the reactions to events.
#' @param evt_name The current event being processed
#' @param input_list_trt A list of simulation inputs
#'
#' @return The modified input list with the updates after executing the corresponding reactions
#'
#' @examples
#' eval_reactevt(x = input_list_trt$evt_react_list,evt_name ="evt1",input_list_trt=input_list_trt)
#'
#' @keywords internal
#' @noRd

eval_reactevt <-  function(x,evt_name,input_list_trt=NULL){
  # Initial set-up --------------------------

  position <- which(evt_name==names(x))
  pos_l <- length(position)
  if (pos_l==0 | pos_l>1 ) {
    stop("Reaction to event ", evt_name, " not recognised or more than one reaction found. Make sure that only one reaction has been defined for the event")    
  }
    
    #Set up values for easier access
    drq <- input_list_trt[["drq"]]
    drc <- input_list_trt[["drc"]]
    prevtime <- input_list_trt[["prevtime"]]
    curtime <- input_list_trt[["curtime"]]
    
    #For each cost/utility category, get the undiscounted/discounted outcomes and get the inputs for ongoing/instant/cycle
    #We do undiscounted first as the discounted value will be overwritten to save some writing time
    
    # Costs -------------------------------------------------------------------
    
    for (cost_cat in input_list_trt$uc_lists$cost_categories_ongoing) {
      input_list_trt[paste0(cost_cat,"_","ongoing_undisc")] <- AddOngoing(lcldr=0,
                                                                          lclprvtime=prevtime,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(cost_cat,"_","ongoing")]])
      input_list_trt[paste0(cost_cat,"_","ongoing")] <- AddOngoing(lcldr=drc,
                                                                   lclprvtime=prevtime,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(cost_cat,"_","ongoing")]])

      input_list_trt[["itemcosts"]] <- input_list_trt[["itemcosts"]] + input_list_trt[[paste0(cost_cat,"_","ongoing")]]
      input_list_trt[["itemcosts_undisc"]] <- input_list_trt[["itemcosts_undisc"]] + input_list_trt[[paste0(cost_cat,"_","ongoing_undisc")]]
      
    }
    
    for (cost_cat in input_list_trt$uc_lists$cost_categories_instant) {
      input_list_trt[paste0(cost_cat,"_","instant_undisc")] <- AddInstant(lcldr=0,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(cost_cat,"_","instant")]])
      input_list_trt[paste0(cost_cat,"_","instant")] <- AddInstant(lcldr=drc,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(cost_cat,"_","instant")]])
      
      input_list_trt[["itemcosts"]] <- input_list_trt[["itemcosts"]] + input_list_trt[[paste0(cost_cat,"_","instant")]]
      input_list_trt[["itemcosts_undisc"]] <- input_list_trt[["itemcosts_undisc"]] + input_list_trt[[paste0(cost_cat,"_","instant_undisc")]]
    }
    
    for (cost_cat in input_list_trt$uc_lists$cost_categories_cycle) {
      if (length(input_list_trt[paste0(cost_cat,"_","cycle")])==1 & input_list_trt[[paste0(cost_cat,"_","cycle")]]==0) {
        input_list_trt[paste0(cost_cat,"_","cycle")] <- list(addcycle=0)
        input_list_trt[paste0(cost_cat,"_","cycle_undisc")] <- list(addcycle=0)
        
        input_list_trt[["itemcosts"]] <- input_list_trt[["itemcosts"]] + input_list_trt[[paste0(cost_cat,"_","cycle")]]
        input_list_trt[["itemcosts_undisc"]] <- input_list_trt[["itemcosts_undisc"]] + input_list_trt[[paste0(cost_cat,"_","cycle_undisc")]]
      } else{
        input_list_trt[paste0(cost_cat,"_","cycle_undisc")] <- AddCycle(lcldr=0,
                                                                        lclprvtime=prevtime,
                                                                        cyclelength = input_list_trt[[paste0(cost_cat,"_","cycle_l")]],
                                                                        lclcurtime=curtime,
                                                                        lclval= input_list_trt[[paste0(cost_cat,"_","cycle")]],
                                                                        starttime = input_list_trt[[paste0(cost_cat,"_","cycle_starttime")]]) #cycles of 1 week
        input_list_trt[paste0(cost_cat,"_","cycle")] <- AddCycle(lcldr=drc,
                                                                 lclprvtime=prevtime,
                                                                 cyclelength = input_list_trt[[paste0(cost_cat,"_","cycle_l")]],
                                                                 lclcurtime=curtime,
                                                                 lclval= input_list_trt[[paste0(cost_cat,"_","cycle")]],
                                                                 starttime = input_list_trt[[paste0(cost_cat,"_","cycle_starttime")]]) #cycles of 1 week
        
        input_list_trt[["itemcosts"]] <- input_list_trt[["itemcosts"]] + input_list_trt[[paste0(cost_cat,"_","cycle")]]
        input_list_trt[["itemcosts_undisc"]] <- input_list_trt[["itemcosts_undisc"]] + input_list_trt[[paste0(cost_cat,"_","cycle_undisc")]]
      }
    }

    # Utilities ----------------------------------------------------------
    
    for (util_cat in input_list_trt$uc_lists$util_categories_ongoing) {
      input_list_trt[paste0(util_cat,"_","ongoing_undisc")] <- AddOngoing(lcldr=0,
                                                                          lclprvtime=prevtime,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(util_cat,"_","ongoing")]])
      input_list_trt[paste0(util_cat,"_","ongoing")] <- AddOngoing(lcldr=drq,
                                                                   lclprvtime=prevtime,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(util_cat,"_","ongoing")]])
      
      input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","ongoing")]]
      input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","ongoing_undisc")]]
      
    }
    
    for (util_cat in input_list_trt$uc_lists$util_categories_instant) {
      input_list_trt[paste0(util_cat,"_","instant_undisc")] <- AddInstant(lcldr=0,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(util_cat,"_","instant")]])
      input_list_trt[paste0(util_cat,"_","instant")] <- AddInstant(lcldr=drq,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(util_cat,"_","instant")]])
      
      input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","instant")]]
      input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","instant_undisc")]]
    }
    
    for (util_cat in input_list_trt$uc_lists$util_categories_cycle) {
      if (length(input_list_trt[[paste0(util_cat,"_","cycle")]])==1 & input_list_trt[[paste0(util_cat,"_","cycle")]]==0) {
        input_list_trt[paste0(util_cat,"_","cycle")] <- list(addcycle=0)
        input_list_trt[paste0(util_cat,"_","cycle_undisc")] <- list(addcycle=0)
        
        input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","cycle")]]
        input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","cycle_undisc")]]
      } else{
        input_list_trt[paste0(util_cat,"_","cycle_undisc")] <- AddCycle(lcldr=0, 
                                                                        lclprvtime=prevtime, 
                                                                        cyclelength = input_list_trt[p[aste0(util_cat,"_","cycle_l")]],
                                                                        lclcurtime=curtime,
                                                                        lclval= input_list_trt[[paste0(util_cat,"_","cycle")]],
                                                                        starttime = input_list_trt[[paste0(util_cat,"_","cycle_starttime")]]) #cycles of 1 week
        input_list_trt[paste0(util_cat,"_","cycle")] <- AddCycle(lcldr=drq,
                                                                 lclprvtime=prevtime,
                                                                 cyclelength = input_list_trt[p[aste0(util_cat,"_","cycle_l")]],
                                                                 lclcurtime=curtime,
                                                                 lclval= input_list_trt[[paste0(util_cat,"_","cycle")]],
                                                                 starttime = input_list_trt[[paste0(util_cat,"_","cycle_starttime")]]) #cycles of 1 week
        
        input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","cycle")]]
        input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","cycle_undisc")]]
      }
    }
    
    additional_ly <- AddOngoing(lcldr=drq,
                                lclprvtime=prevtime, 
                                lclcurtime=curtime,
                                lclval=1)
    additional_ly_undisc <- AddOngoing(lcldr=0,
                                       lclprvtime=prevtime,
                                       lclcurtime=curtime,
                                       lclval=1)
    
    input_list_trt[["itemlys"]] <- additional_ly
    input_list_trt[["itemlys_undisc"]] <- additional_ly_undisc
    
    
    #evaluate reaction
    
    input_list_trt <- local({
      input_list_trt <- input_list_trt
      eval(x[[position]][["react"]]) #run script
      out <- input_list_trt
    },input_list_trt)
    
    
    
  
  return(input_list_trt)
  

}


# Get Inputs for Event -----------------------------------------------------------------------------------------------------------------------------------------------

#' Evaluates the cost/utility/cycle unevaluated expressions to be processed by the simulation engine
#'
#' @param x The specific cost/utility and its type (ongoing, instant...) to be used, created through add_cost/add_util
#' @param ifnull Value to be used if the input has not been defined
#' @param type Identifies what type of input is being used. Can be "cost", "util", "cycle_l" (cycle length) and "cycle_starttime" (starting time of the cycle)
#' @param evt_trt_i The event-intervention identifier to understand which specific input to use, separated by an underscore
#' @param input_list_trt_i  A list of simulation inputs
#'
#' @return A numeric vector of evaluated costs/utilities/cycle lengths/starting times for the specific event and intervention defined
#'
#' @examples
#' get_input(x = input_list_trt$uc_lists$cost_ongoing_list,ifnull=0,type="cost",evt_trt_i="evt1_int",input_list_trt_i=input_list_trt)
#'
#' @keywords internal
#' @noRd

get_input <-  function(x,ifnull=0,type,evt_trt_i =evt_trt, input_list_trt_i=input_list_trt){
  out <- NULL
  items <- x[names(x)==evt_trt_i]
  items_l <- length(items)
  if (items_l==0) {
    out <-  c(out,ifnull)
  } else{
    #If length is 1, then don't do loop as it's slow
    if (items_l==1) {
      out <-  c(out,
                if(is.null(items[[1]][[type]])){
                  ifnull
                } else{
                  eval(items[[1]][[type]],input_list_trt_i)
                } #lazy eval will give error on null, so just put 0 in that case
      )
      #otherwise do loop per item
    } else{
      for (i in 1:items_l) {
        out <-  c(out,
                  if(is.null(items[[i]][[type]])){
                    ifnull
                  } else{
                    eval(items[[i]][[type]],input_list_trt_i)
                  } #lazy eval will give error on null, so just put 0 in that case
        )
      }
    }

  }
  return(out)

}


# Helper function to create mean and 95% interval -------------------------

#' Calculate mean and 95% CI from samples
#'
#' @param x  The output_psa data frame from the list object returned by `RunSim()`
#' @param element Variable for which mean and 95% CIs are computed (single string)
#' @param trt Treatment for which mean and 95% CIs are computed (single string)
#' @param round_digit Number of digits to round outputs
#'
#' @return Mean and 95% CI from the PSA samples
#' @importFrom purrr map_dbl
#' @importFrom stats quantile
#'
#' @examples
#' interval_out(x=results$output_psa,element="costs.",trt="int",round_digit=3)
#'
#' @keywords internal
#' @noRd

interval_out <- function(x, element, trt,round_digit=2) {
  out <- paste0(round(mean(map_dbl(x,paste0(element,trt)),na.rm=TRUE),round_digit),
                "(",
                round(quantile(map_dbl(x,paste0(element,trt)),0.025,na.rm=TRUE),round_digit) ,
                ", ",
                round(quantile(map_dbl(x,paste0(element,trt)),0.975,na.rm=TRUE),round_digit),
                ")"
  )
  return(out)
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
#' @keywords internal
#' @noRd

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
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
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
