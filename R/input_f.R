# Replicate profiles --------------------------------------------------------
#' Replicate profiles data.frame
#'
#' @param profiles data.frame of profiles
#' @param replications integer, final number of observations
#' @param probabilities vector of probabilities with the same length as the number of rows of profiles. Does not need to add up to 1 (are reweighted)
#' @param replacement Boolean whether replacement is used
#' @param seed_used Integer with the seed to be used for consistent results
#'
#' @return Resampled data.frame of profiles
#' @export
#' 
#'
#' @examples
#' \dontrun{
#' replicate_profiles(profiles=data.frame(id=1:100,age=rnorm(100,60,5)),
#' replications=200,probabilities=rep(1,100))
#' }
replicate_profiles <- function(profiles,
                               replications,
                               probabilities = NULL,
                               replacement = TRUE,
                               seed_used = NULL
){
  
  if (!is.null(seed_used)) {
    set.seed(seed_used)
  }
  
  rows_sampled <- sample.int(nrow(profiles),size=replications,prob=probabilities,replace=replacement)
  output <- profiles[rows_sampled,]
  rownames(output) <- NULL  
  
  return(output)
}

# Create indicators for sensitivity/DSA analysis --------------------------------------------------------
#' Creates a vector of indicators (0 and 1) for sensitivity/DSA analysis
#'
#' @param sens current analysis iterator
#' @param n_sensitivity total number of analyses to be run
#' @param elem vector of 0s and 1s of elements to iterate through (1 = parameter is to be included in scenario/DSA)
#' @param n_elem_before Sum of 1s (# of parameters to be included in scenario/DSA) that go before elem
#'
#' @return Numeric vector composed of 0 and 1, where value 1 will be used by `pick_val_v` to pick the corresponding index in its `sens` argument
#' @export
#' 
#' @details
#' n_elem_before is to be used when several indicators want to be used (e.g., for patient level and common level inputs) while facilitating readibility of the code
#'
#' @examples
#' \dontrun{
#'create_indicators(10,20,c(1,1,1,1))
#'create_indicators(7,20,c(1,0,0,1,1,1,0,0,1,1),2)
#' }
create_indicators <- function(sens,n_sensitivity,elem,n_elem_before=0){
  
  if (n_sensitivity<sens) {
    stop("n_sensitivity is smaller than the iterator sens")
  }
  
  if (any(elem>1| elem<0)) {
    stop("One or more elements of elem is >1 or <0. It should be a vector of 0s and 1s.")
  }
  n_elem <- sum(elem)

    if (n_elem_before + n_elem >n_sensitivity) {
    stop("n_sensitivity is smaller than n_elem_before + sum(elem)")
  }

  n_elem_total <- length(elem)
  
  if ((n_elem + n_elem_before <sens) | (n_elem_before >=sens)) { #here fix
    out <- rep(0, n_elem_total)
  } else{
    # which position to use to put the value 1 in indicator
    pos_indicator <-  sens - n_sensitivity*floor((sens-1)/n_sensitivity)  - n_elem_before
    pos_indicator <- which(elem==1)[pos_indicator]
    
    out <- append(rep(0, n_elem_total)[-pos_indicator],1,pos_indicator-1) 
  }
  return(out)
}


# Helper to draw PSA values --------------------------------------------------------
#' Helper function to create a list with random draws or whenever a series of functions needs to be called. Can be implemented within `pick_val_v`.
#'
#' @param f A string or vector of strings with the function to be called, e.g., "rnorm"
#' @param ... parameters to be passed to the function (e.g., if "rnorm", arguments `n`, `mean`, `sd`)
#'
#' @return List with length equal to `f` of parameters called
#' @export
#' 
#' @details
#' This function can be used to pick values for the PSA within `pick_val_v.` 
#' 
#' The function will ignore NA items within the respective parameter (see example below).
#' If an element in f is NA (e.g., a non PSA input) then it will return NA as its value
#' This feature is convenient when mixing distributions with different number of arguments, e.g., `rnorm` and `rgengamma`.
#' 
#' While it's slightly lower than individually calling each function, it makes the code easier to read and more transparent
#'
#' @examples
#' \dontrun{
#' params <- list(
#' param=list("a","b"),
#' dist=list("rlnorm","rnorm"),
#' n=list(4,1),
#' a=list(c(1,2,3,4),1),
#' b=list(c(0.5,0.5,0.5,0.5),0.5),
#' dsa_min=list(c(1,2,3,4),2),
#' dsa_max=list(c(1,2,3,4),3)
#' )
#'pick_psa(params[["dist"]],params[["n"]],params[["a"]],params[["b"]])
#'
#'#It works with functions that require different number of parameters
#'params <- list(
#'  param=list("a","b","c"),
#'  dist=list("rlnorm","rnorm","rgengamma"),
#'  n=list(4,1,1),
#'  a=list(c(1,2,3,4),1,0),
#'  b=list(c(0.5,0.5,0.5,0.5),0.5,1),
#'  c=list(NA,NA,0.2),
#'  dsa_min=list(c(1,2,3,4),2,1),
#'  dsa_max=list(c(1,2,3,4),3,3)
#')
#'
#'pick_psa(params[["dist"]],params[["n"]],params[["a"]],params[["b"]],params[["c"]])
#'
#' #Can be combined with multiple type of functions and distributions if parameters are well located
# params <- list(
# param=list("a","b","c","d"),
# dist=list("rlnorm","rnorm","rgengamma","draw_tte"),
# n=list(4,1,1,1),
# a=list(c(1,2,3,4),1,0,"norm"),
# b=list(c(0.5,0.5,0.5,0.5),0.5,1,1),
# c=list(NA,NA,0.2,0.5), #NA arguments will be ignored
# c=list(NA,NA,NA,NA), #NA arguments will be ignored
# dsa_min=list(c(1,2,3,4),2,1,0),
# dsa_max=list(c(1,2,3,4),3,3,2)
# )
#'
#'#'params <- list(
#' param=list("a","b","c","d"),
#' dist=list("rlnorm","rnorm","rgengamma","draw_tte"),
#' n=list(4,1,1,1),
#' a=list(c(1,2,3,4),1,0,"norm"),
#' b=list(c(0.5,0.5,0.5,0.5),0.5,1,1),
#' c=list(NA,NA,0.2,0.5),
#' c=list(NA,NA,NA,NA), #NA arguments will be ignored
#' dsa_min=list(c(1,2,3,4),2,1,0),
#' dsa_max=list(c(1,2,3,4),3,3,2)
#' )
#'
#'
#'
#' }
#' 
pick_psa <- function(f,...){
  args_called <- list(...)
  sapply(1:length(f), function(x) {
    if(!is.na(f[[x]])){
    args <- lapply(args_called, `[[`, x)
    args <- args[!sapply(args, function(x) any(is.na(x)))]
    do.call(f[[x]], args)
    } else{
      NA
    }
  })
}

# Select which values to apply --------------------------------------------------------
#' Select which values should be applied in the corresponding loop for several values (vector or list).
#'
#' @param base Value if no PSA/DSA/Scenario
#' @param psa Value if PSA
#' @param sens Value if DSA/Scenario
#' @param psa_ind Boolean whether PSA is active
#' @param sens_ind Boolean whether Scenario/DSA is active
#' @param indicator Indicator which checks whether the specific parameter/parameters is/are active in the DSA or Scenario loop 
#' @param indicator_psa Indicator which checks whether the specific parameter/parameters is/are active in the PSA loop.
#'  If NULL, it's assumed to be a vector of 1s of length equal to length(indicator)
#' @param names_out Names to give the output list
#'
#' @return List of used for the inputs
#' @export
#' 
#' @details
#' This function can be used with vectors or lists, but will always return a list.
#' Lists should be used when correlated variables are introduced to make sure the selector knows how to choose among those
#'
#' @examples
#' \dontrun{
#' pick_val_v(base = c(0,0),
#'            psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
#'            sens = c(2,3),
#'            psa_ind = TRUE,
#'            sens_ind = FALSE,
#'            indicator=c(1,0)
#'            )
#'            
#' pick_val_v(base = c(0,0),
#' psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
#' sens = c(2,3),
#' psa_ind = TRUE,
#' sens_ind = TRUE,
#' indicator=c(1,0),
#' indicator_psa=c(0,1)
#' )
#' pick_val_v(base = c(c(2,3),list(c(1,2))),
#'             psa =sapply(1:3,
#'                         function(x) eval(call(
#'                           c("rnorm","rnorm","mvrnorm")[[x]],
#'                           1,
#'                           c(c(2,3),list(c(1,2)))[[x]],
#'                           c(c(0.1,0.1),list(matrix(c(1,0.1,0.1,1),2,2)))[[x]]
#'                         ))),
#'             sens = c(c(2,3),list(c(1,2))),
#'             psa_ind = TRUE,
#'             sens_ind = FALSE,
#'             indicator=c(1,0,0),
#'             names_out=c("util","util2","correlated_vector")
#'            )
#' 
#' pick_val_v(base        = df_par[,"base_value"],
#'            psa         = sapply(1:nrow(df_par), function(x)
#'                            eval(call(df_par[x,"PSA_dist"],1,df_par[x,"a"],df_par[x,"b"]))),
#'            sens        = df_par[,sensitivity_names[sens_sel]],
#'            psa_ind     = TRUE, sens_ind = FALSE,
#'            indicator   = indicators,
#'            names_out   = df_par[,"parameter_name"]
#'            )
#'            
#'  pick_val_v(
#'   base = c(0,0),
#'   psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
#'   sens = c(2,3),
#'   psa_ind = TRUE,
#'   sens_ind = TRUE,
#'   indicator=create_indicators(sens=1,n_sensitivity=2,n_elem = 2)
#'  )
#'  
#' }
pick_val_v <- function(base,
                       psa,
                       sens,
                       psa_ind = psa_bool,
                       sens_ind = sens_bool,
                       indicator,
                       indicator_psa = NULL,
                       names_out=NULL
){


  if ((any(!indicator %in% c(0,1))) | (!psa_ind %in% c(0,1)) | (!sens_ind %in% c(0,1)) ) {
    stop("Indicator, psa_ind or sens_ind are not FALSE/TRUE (or 0/1)")
  }
  
  len_ind <- length(indicator)
  
  if (is.null(indicator_psa)) {
    indicator_psa <- rep(1,len_ind)
  }else{
    if(len_ind != length(indicator_psa)){
      stop("Length of indicator vector is different than length of indicator_psa")
    }
  }
  
  output <- list()
    #if the parameter is out of the dsa/scenario specific iteration, or not in DSA/scenario, use PSA/normal. 
  for (it in 1:len_ind) {
    output[[it]] <-  if (indicator[it]==0 | sens_ind==F ) { 
      if (psa_ind==T & indicator_psa[it]==1) {
        psa[[it]]
        } else {
          base[[it]]
          }
      } else {#If active, use DSA if DSA and scenario if scenario
        sens[[it]]
      } 
  }
  
  if (!is.null(names_out)) {
    names(output) <- names_out
  }
    
  return(output)
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
#' Whenever a function is directly implemented which must be evaluated later and that has no object name attached (e.g., `pick_val_v`), 
#' it should be implemented after a first `add_item()` (empty or with content) to avoid confusing the `.data` argument, or wrapping the function within `substitute()`
#'
#' @examples
#' \dontrun{
#' add_item(fl.idfs = 0)
#' add_item(util_idfs = if(psa_bool){rnorm(1,0.8,0.2)} else{0.8}, util.mbc = 0.6, cost_idfs = 2500)
#' common_inputs <- add_item() %>%
#' add_item(pick_val_v(
#'   base      = l_statics[["base"]],
#'   psa       = pick_psa(
#'     l_statics[["function"]],
#'     l_statics[["n"]],
#'     l_statics[["a"]],
#'     l_statics[["b"]]
#'   ),
#'   sens      = l_statics[[sens_name_used]],
#'   psa_ind   = psa_bool,
#'   sens_ind  = sensitivity_bool,
#'   indicator = indicators_statics,
#'   names_out = l_statics[["parameter_name"]]
#' )
#' )
#' }
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

new_event <- function(evt){
  new_evt_name <- names(evt)
  new_evt <- setNames(unlist(evt),new_evt_name)
  if (!is.numeric(new_evt)) {
    stop("New event times are not all numeric, please review")
  }
  
  input_list_arm <- parent.frame()$input_list_arm
  
  evtlist_temp <- list(cur_evtlist = c(input_list_arm$cur_evtlist,
                                  new_evt))
  
  input_list_arm[["cur_evtlist"]] <- evtlist_temp$cur_evtlist
  
  list2env(input_list_arm["cur_evtlist"],envir = parent.frame())
  assign("input_list_arm",input_list_arm, envir = parent.frame())

}

# Modify event in list of events ---------------------------------------------

#' Modify the time of existing events
#'
#' @param evt A list of events and their times
#' @param create_if_null A boolean.
#'  If TRUE, it will create non-existing events with the chosen time to event.
#'  If FALSE, it will ignore those.
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

modify_event <- function(evt,create_if_null=TRUE){
  input_list_arm <- parent.frame()$input_list_arm
  
  evt_unlist <- unlist(evt)
  if(!is.numeric(evt_unlist)){
      stop("Modify event times are not all numeric, please review")
  }
  names_evt <- names(evt)
  
  if (create_if_null==FALSE) {
    names_obj_temp <- names(input_list_arm$cur_evtlist)
    names_found <- names_evt[names_evt %in% names_obj_temp]
    if (length(names_found)==0) {
      warning("Some or all event/s in modify_evt within ", paste(names_evt,collapse=", "), " not found. Use new_evt to add new events.")
    }
    matched <- which(names_obj_temp %in% names_found)
    
    input_list_arm[["cur_evtlist"]][matched] <- evt_unlist[names_obj_temp[names_obj_temp %in% names_found]]
  } else{
    input_list_arm[["cur_evtlist"]][names_evt] <- evt_unlist
  }
  
   
  list2env(input_list_arm["cur_evtlist"],envir = parent.frame())
  assign("input_list_arm",input_list_arm, envir = parent.frame())
  

}





# Modify item in input list -------------------------------------------------------------------------------------------------------------------------------

#' Modify the value of existing items
#'
#' @param list_item A list of items and their values or expressions
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_item` with a list of one element, it's better to group them into a single `modify_item` with a list of two elements.
#'
#' Costs and utilities can be modified by using the construction `type_name_category`, where type is either "qaly" or "cost",
#'  name is the name (e.g., "default") and category is the category used (e.g., "instant"), so one could pass `cost_default_instant` and modify the cost.
#'  This will overwrite the value defined in the corresponding cost/utility section.
#'
#' @examples
#' \dontrun{
#' modify_item(list(cost.idfs = 500, cost.tx = cost.tx + 4000))
#' }

modify_item <- function(list_item){
  input_list_arm <- parent.frame()$input_list_arm
  
  input_list_arm[names(list_item)] <- lapply(list_item, unname)
  
  list2env(list_item,envir = parent.frame())
  assign("input_list_arm",input_list_arm, envir = parent.frame())
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
#' \dontrun{
#' add_reactevt(name_evt = "start",input = {})
#' add_reactevt(name_evt = "idfs",input = {modify_item(list("fl.idfs"= 0))})
#' }

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
#' @param arm The intervention for which the events and initial event times are defined
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
#' \dontrun{
#' add_tte(arm="int",evts = c("start","ttot","idfs","os"),
#' input={
#' start <- 0
#' idfs <- draw_tte(1,'lnorm',coef1=2, coef2=0.5)
#' ttot <- min(draw_tte(1,'lnorm',coef1=1, coef2=4),idfs)
#' os <- draw_tte(1,'lnorm',coef1=0.8, coef2=0.2)
#' })
#' }
#'
add_tte <- function(.data=NULL,arm, evts, other_inp = NULL,input){
  data_list <- .data

  for (arm_i in arm) {

    if (!is.character(other_inp) & !is.null(other_inp)) {
      stop("other_inp argument is required to be a character vector or be set to NULL")
    }

    if (!is.character(evts) | length(evts)<2) {
      stop("evts argument in add_tte for the intervention ", arm_i, " is required to be a character vector with length >1")
    }


    evt_l <- list(list(expr=substitute(input),
                       evts = evts,
                       other_inp = other_inp
    ))
    names(evt_l) <- paste(arm_i)

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
#' @examples \dontrun{
#' disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' }
#'
#' @export

disc_ongoing <- function(lcldr=0.035, lclprvtime, lclcurtime, lclval){
  
  Instantdr <- log(1+lcldr)
  
  # calculate additional qalys
  if (is.null(lclval)) {
    add <- 0
  } else if(lcldr==0) {
    add <- lclval*(lclcurtime - lclprvtime)
  } else{
    add <- ((lclval)/(0 - Instantdr)) * (exp(lclcurtime * ( 0 - Instantdr)) - exp(lclprvtime * (0 - Instantdr)))
    
  }
  
  return(add)
}


#' Calculate discounted costs and qalys between events for vectors
#'
#' @param lcldr The discount rate 
#' @param lclprvtime The time of the previous event in the simulation
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on continuous time discounting
#'
#' @examples \dontrun{
#' disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' }
#'
#' @export

disc_ongoing_v <- function(lcldr=0.035, lclprvtime, lclcurtime, lclval){
  
  Instantdr <- log(1+lcldr)
  
  # calculate additional qalys
  if(lcldr==0) {
    add <- lclval*(lclcurtime - lclprvtime)
  } else{
    add <- ((lclval)/(0 - Instantdr)) * (exp(lclcurtime * ( 0 - Instantdr)) - exp(lclprvtime * (0 - Instantdr)))
    
  }
  
  return(add)
}

#' Calculate instantaneous discounted costs or qalys for vectors
#'
#' @param lcldr The discount rate
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on discrete time discounting
#'
#' @examples \dontrun{
#' disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500)
#' }
#'
#' @export
#' 
disc_instant_v <- function(lcldr=0.035, lclcurtime, lclval){
  

    addinst <- lclval * ((1+lcldr)^(-lclcurtime))    

  return(addinst)
}

#' Calculate instantaneous discounted costs or qalys
#'
#' @param lcldr The discount rate
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on discrete time discounting
#'
#' @examples \dontrun{
#' disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500)
#' }
#'
#' @export
#' 
disc_instant <- function(lcldr=0.035, lclcurtime, lclval){
  
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
#' @examples \dontrun{
#' disc_cycle(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' }
#'
#' @export

disc_cycle <- function(lcldr=0.035, lclprvtime=0, cyclelength,lclcurtime, lclval,starttime=0){
  
  addcycle <- 0
  
  if (is.null(lclval) ) {
    addcycle <- 0
  } else{
    
    #Note this makes the cycle utilities work weird, so do not use cycle utilities for now!
    for (i in 1:length(lclval)) {
      lclval_i <- lclval[i]
      starttime_i <- starttime[i]
      cyclelength_i <- cyclelength[i]
      
      if (lclval_i==0 ) {} else{
        
        cycle.time.total <- if(starttime_i>= lclcurtime){0}else{seq(from=starttime_i, to = lclcurtime , by= cyclelength_i)} #all cycles that happened until current time of event
        
        #If the cost starts at the selected starttime or at time 0, then include that time, otherwise exclude it
        if (lclprvtime==0) {
          cycle.time <- c(0,cycle.time.total[cycle.time.total >= lclprvtime])  #times at which the cycles take place during this event, put this condition to count also time 0
          n_cycles <- length(unique(cycle.time))
          s <- (1+lcldr)^cyclelength_i -1
          if (lcldr==0) {
            addcycle <- sum(addcycle, lclval_i * n_cycles)
          } else{
            addcycle <- sum(addcycle,lclval_i * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1) )
          }
        } else{
          if (starttime_i ==lclprvtime & lclprvtime==0) {
            cycle.time <- cycle.time.total[cycle.time.total >= lclprvtime]  #times at which the cycles take place during this event, put this condition to count also time of the previous event
          } else{
            cycle.time <- cycle.time.total[cycle.time.total > lclprvtime]  #times at which the cycles take place during this event
          }
          n_cycles_remaining <- length(cycle.time)
          d <- lclprvtime/cyclelength_i
          s <- (1+lcldr)^cyclelength_i -1
          if (lcldr==0) {
            addcycle <- sum(addcycle, lclval_i * n_cycles_remaining)
          } else{
            addcycle <- sum(addcycle, lclval_i * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d)) )
          }
          
        }
        
        #If starting from 0, can be changed substituting interest rate such that s = (1+r)^cyclelength - 1, and using the formula that lclvalq * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1)
        #If starting from time t, then compute transformed time as d = t/cyclelength and use lclvalq * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d-1)), where
        #n_cycles_remaining is the n_cycles - d (so the remaining cycles to be considered), e.g. if 13 cycles (From t=0), and delay 6 periods, then n_cycles_remaining = 7 and d=6
        
      }
      
    }
  }
  return(addcycle)
}

#' Cycle discounting for vectors
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
#' @examples \dontrun{
#' disc_cycle_v(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' }
#'
#' @export

disc_cycle_v <- function(lcldr=0.035, lclprvtime=0, cyclelength,lclcurtime, lclval,starttime=0){
  
  addcycle <- rep(0,length(lclval))
  
  #Note this makes the cycle utilities work weird, so do not use cycle utilities for now!
  for (i in 1:length(lclval)) {
    lclval_i <- lclval[i]
    starttime_i <- starttime[i]
    cyclelength_i <- cyclelength[i]
    lclcurtime_i <- lclcurtime[i]
    lclprvtime_i <- lclprvtime[i]
    
    
    if (lclval_i==0 ) {} else{
      
      cycle.time.total <- if(starttime_i>= lclcurtime_i){0}else{seq(from=starttime_i, to = lclcurtime_i , by= cyclelength_i)} #all cycles that happened until current time of event
      
      #If the cost starts at the selected starttime or at time 0, then include that time, otherwise exclude it
      if (lclprvtime_i==0) {
        cycle.time <- c(0,cycle.time.total[cycle.time.total >= lclprvtime_i])  #times at which the cycles take place during this event, put this condition to count also time 0
        n_cycles <- length(unique(cycle.time))
        s <- (1+lcldr)^cyclelength_i -1
        if (lcldr==0) {
          addcycle[i] <- sum(addcycle[i],lclval_i * n_cycles )
        } else{
          addcycle[i] <- sum(addcycle[i],lclval_i * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1) )
        }
      } else{
        if (starttime_i ==lclprvtime_i & lclprvtime_i==0) {
          cycle.time <- cycle.time.total[cycle.time.total >= lclprvtime_i]  #times at which the cycles take place during this event, put this condition to count also time of the previous event
        } else{
          cycle.time <- cycle.time.total[cycle.time.total > lclprvtime_i]  #times at which the cycles take place during this event
        }
        n_cycles_remaining <- length(cycle.time)
        d <- lclprvtime_i/cyclelength_i
        s <- (1+lcldr)^cyclelength_i -1
        if (lcldr==0) {
          addcycle[i] <- sum(addcycle[i], lclval_i * n_cycles_remaining )
        } else{
          addcycle[i] <- sum(addcycle[i], lclval_i * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d)) )
        }
      }
      
      #If starting from 0, can be changed substituting interest rate such that s = (1+r)^cyclelength - 1, and using the formula that lclvalq * (1 - (1+s)^-n_cycles)/(s*(1+s)^-1)
      #If starting from time t, then compute transformed time as d = t/cyclelength and use lclvalq * (1 - (1+s)^-n_cycles_remaining)/(s*(1+s)^(d-1)), where
      #n_cycles_remaining is the n_cycles - d (so the remaining cycles to be considered), e.g. if 13 cycles (From t=0), and delay 6 periods, then n_cycles_remaining = 7 and d=6
    }
    
    
  }
  return(addcycle)
}


#' Reverts list
#'
#' @param ls List
#'
#' @return Reverted list
#'
#' @examples \dontrun{
#' revert_list(ls=list(a=list(ab=1,ac=2),b=list(ab=2,ac=3)))
#' }
#'
#' @export

revert_list <- function(ls) {
  if(!is.list(ls)){stop("Object passed is not a list")}
  
  x <- lapply(ls,`[`, names(ls[[1]]))
  apply(do.call(rbind, x), 2, as.list)
}
