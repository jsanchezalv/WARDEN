# Initial event list --------------------------------------------------------------------------------------------------------------------------------------

#' Execute the initial time to events and separate the events from other inputs that are stored
#'
#' @param trt_name A character string of the name of the intervention
#' @param input_list_trt A list of simulation inputs
#'
#' @return A named vector of initial event times, and a named vector of other inputs to be stored
#'
#' @examples
#' initiate_evt(trt = "int",input_list_trt = input_list_trt)
#'
#' @keywords internal
#' @noRd


initiate_evt <- function(trt_name,input_list_trt){
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
#' get_next_evt(evt_list = input_list_trt$cur_evtlist)
#'
#' @keywords internal
#' @noRd

get_next_evt <- function(evt_list){                  # This function identifies which event is to be processed next for each patient, depending on intervention

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
#' react_evt(thisevt="evt1",trt="int",input_list_trt=input_list_trt)
#'
#' @keywords internal
#' @noRd

react_evt <- function(thisevt,trt,input_list_trt=NULL){      # This function processes the next event (as identified in the GetNextEvt function)
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
                                                                        input_list_trt_i=input_list_trt)
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
                                                                        input_list_trt_i=input_list_trt)
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
      input_list_trt[paste0(cost_cat,"_","ongoing_undisc")] <- disc_ongoing(lcldr=0,
                                                                          lclprvtime=prevtime,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(cost_cat,"_","ongoing")]])
      
      input_list_trt[paste0(cost_cat,"_","ongoing")] <- disc_ongoing(lcldr=drc,
                                                                   lclprvtime=prevtime,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(cost_cat,"_","ongoing")]])

      input_list_trt[["itemcosts"]] <- input_list_trt[["itemcosts"]] + input_list_trt[[paste0(cost_cat,"_","ongoing")]]
      
      input_list_trt[["itemcosts_undisc"]] <- input_list_trt[["itemcosts_undisc"]] + input_list_trt[[paste0(cost_cat,"_","ongoing_undisc")]]
      
    }
    
    for (cost_cat in input_list_trt$uc_lists$cost_categories_instant) {
      input_list_trt[paste0(cost_cat,"_","instant_undisc")] <- disc_instant(lcldr=0,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(cost_cat,"_","instant")]])
      
      input_list_trt[paste0(cost_cat,"_","instant")] <- disc_instant(lcldr=drc,
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
        input_list_trt[paste0(cost_cat,"_","cycle_undisc")] <- disc_cycle(lcldr=0,
                                                                        lclprvtime=prevtime,
                                                                        cyclelength = input_list_trt[[paste0(cost_cat,"_","cycle_l")]],
                                                                        lclcurtime=curtime,
                                                                        lclval= input_list_trt[[paste0(cost_cat,"_","cycle")]],
                                                                        starttime = input_list_trt[[paste0(cost_cat,"_","cycle_starttime")]]) #cycles of 1 week
        
        input_list_trt[paste0(cost_cat,"_","cycle")] <- disc_cycle(lcldr=drc,
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
      input_list_trt[paste0(util_cat,"_","ongoing_undisc")] <- disc_ongoing(lcldr=0,
                                                                          lclprvtime=prevtime,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(util_cat,"_","ongoing")]])
      
      input_list_trt[paste0(util_cat,"_","ongoing")] <- disc_ongoing(lcldr=drq,
                                                                   lclprvtime=prevtime,
                                                                   lclcurtime=curtime,
                                                                   lclval=input_list_trt[[paste0(util_cat,"_","ongoing")]])
      
      input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","ongoing")]]
      
      input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","ongoing_undisc")]]
      
    }
    
    for (util_cat in input_list_trt$uc_lists$util_categories_instant) {
      input_list_trt[paste0(util_cat,"_","instant_undisc")] <- disc_instant(lcldr=0,
                                                                          lclcurtime=curtime,
                                                                          lclval=input_list_trt[[paste0(util_cat,"_","instant")]])
      
      input_list_trt[paste0(util_cat,"_","instant")] <- disc_instant(lcldr=drq,
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
        input_list_trt[paste0(util_cat,"_","cycle_undisc")] <- disc_cycle(lcldr=0, 
                                                                        lclprvtime=prevtime, 
                                                                        cyclelength = input_list_trt[p[aste0(util_cat,"_","cycle_l")]],
                                                                        lclcurtime=curtime,
                                                                        lclval= input_list_trt[[paste0(util_cat,"_","cycle")]],
                                                                        starttime = input_list_trt[[paste0(util_cat,"_","cycle_starttime")]]) #cycles of 1 week
        
        input_list_trt[paste0(util_cat,"_","cycle")] <- disc_cycle(lcldr=drq,
                                                                 lclprvtime=prevtime,
                                                                 cyclelength = input_list_trt[p[aste0(util_cat,"_","cycle_l")]],
                                                                 lclcurtime=curtime,
                                                                 lclval= input_list_trt[[paste0(util_cat,"_","cycle")]],
                                                                 starttime = input_list_trt[[paste0(util_cat,"_","cycle_starttime")]]) #cycles of 1 week
        
        input_list_trt[["itemqalys"]] <- input_list_trt[["itemqalys"]] + input_list_trt[[paste0(util_cat,"_","cycle")]]
        
        input_list_trt[["itemqalys_undisc"]] <- input_list_trt[["itemqalys_undisc"]] + input_list_trt[[paste0(util_cat,"_","cycle_undisc")]]
      }
    }
    
    #LYs
    additional_ly <- disc_ongoing(lcldr=drq,
                                lclprvtime=prevtime, 
                                lclcurtime=curtime,
                                lclval=1)
    
    additional_ly_undisc <- disc_ongoing(lcldr=0,
                                       lclprvtime=prevtime,
                                       lclcurtime=curtime,
                                       lclval=1)
    
    input_list_trt[["itemlys"]] <- additional_ly
    
    input_list_trt[["itemlys_undisc"]] <- additional_ly_undisc
    
    

# Evaluate reaction -------------------------------------------------------

    
    
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
#' @param x  The output_psa data frame from the list object returned by `run_sim()`
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


