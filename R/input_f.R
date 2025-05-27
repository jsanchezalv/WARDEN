
# Global variables for CRAN check -----------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      c('psa_bool',
        'sens_bool',
        'evt_arm',
        'input_list_arm',
        'new_event_name')
    )) 
}

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
#' replicate_profiles(profiles=data.frame(id=1:100,age=rnorm(100,60,5)),
#' replications=200,probabilities=rep(1,100))
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
#'create_indicators(10,20,c(1,1,1,1))
#'create_indicators(7,20,c(1,0,0,1,1,1,0,0,1,1),2)

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


# Create an iterator based on sens of the current iteration within a scenario (DSA) --------------------------------------------------------
#' Create an iterator based on sens of the current iteration within a scenario (DSA)
#'
#' @param sens current analysis iterator
#' @param n_sensitivity total number of analyses to be run
#'
#' @return Integer iterator based on the number of sensitivity analyses being run and the total iterator
#' @export
#' 
#' @details
#' In a situation like a DSA, where two (low and high) scenarios are run, sens will go from 1 to n_sensitivity*2. However,
#' this is not ideal as the parameter selector may depend on knowing the parameter order (i.e., 1, 2, 3...), which means
#' resetting the counter back to 1 once sens reaches n_sensitivity (or any multiple of n_sensitivity) is needed.
#'
#' @examples
#' sens_iterator(5,20)
#' sens_iterator(25,20)

sens_iterator <- function(sens,n_sensitivity){
  
  out <- sens - n_sensitivity*floor((sens-1)/n_sensitivity)
  
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
#' params <- list(
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
#' @param indicator_sens_binary Boolean, TRUE if parameters will be varied fully, FALSE if some elements of the parameters may be changed but not all
#' @param sens_iterator Current iterator number of the DSA/scenario being run, e.g., 5 if it corresponds to the 5th DSA parameter being changed
#' @param distributions List with length equal to length of base where the distributions are stored
#' @param covariances List with length equal to length of base where the variance/covariances are stored (only relevant if multivariate normal are being used)
#'
#' @return List used for the inputs
#' @export
#' 
#' @details
#' This function can be used with vectors or lists, but will always return a list.
#' Lists should be used when correlated variables are introduced to make sure the selector knows how to choose among those
#' This function allows to choose between using an approach where only the full parameters are varied, and an approach where subelements of the parameters can be changed
#'
#' @examples
#' pick_val_v(base = list(0,0),
#'              psa =list(rnorm(1,0,0.1),rnorm(1,0,0.1)),
#'              sens = list(2,3),
#'              psa_ind = FALSE,
#'              sens_ind = TRUE,
#'              indicator=list(1,2),
#'              indicator_sens_binary = FALSE,
#'              sens_iterator = 2,
#'              distributions = list("rnorm","rnorm")
#' )
#' 
#' pick_val_v(base = list(2,3,c(1,2)),
#'              psa =sapply(1:3,
#'                          function(x) eval(call(
#'                            c("rnorm","rnorm","mvrnorm")[[x]],
#'                            1,
#'                            c(2,3,list(c(1,2)))[[x]],
#'                            c(0.1,0.1,list(matrix(c(1,0.1,0.1,1),2,2)))[[x]]
#'                          ))),
#'              sens = list(4,5,c(1.3,2.3)),
#'              psa_ind = FALSE,
#'              sens_ind = TRUE,
#'              indicator=list(1,2,c(3,4)),
#'              names_out=c("util","util2","correlated_vector") ,
#'              indicator_sens_binary = FALSE,
#'              sens_iterator = 4,
#'              distributions = list("rnorm","rnorm","mvrnorm"),
#'              covariances = list(0.1,0.1,matrix(c(1,0.1,0.1,1),2,2))
#' )
#'  
pick_val_v <- function(base,
                         psa,
                         sens,
                         psa_ind = psa_bool,
                         sens_ind = sens_bool,
                         indicator,
                         indicator_psa = NULL,
                         names_out=NULL,
                         indicator_sens_binary = TRUE,
                         sens_iterator = NULL,
                         distributions = NULL,
                         covariances = NULL
){
  
  output <- list()
  
  if(indicator_sens_binary == FALSE & sens_ind == TRUE){ #if the indicators for the sensitivity analyses are not 0/1 but instead a list of integers (elements of parameters)
    
    if(is.null(sens_iterator)|is.null(distributions)){
      stop("sens_iterator nor distributions arguments cannot be NULL if indicator_sens_binary argument is FALSE. sens_indicator should take the value of the corresponding DSA/scenario iterator")
    }
    
    if(psa_ind){
      temp_data <- psa #preassign base value
    }else{
      temp_data <- base #preassign base value
    }
    output <- temp_data #preassign base value
    
    
    #Iterate over each parameter, check which elements are active in the current iterator, then check if an adjustment needs to be done 
    for (i in 1:length(base)) {
      ind <- indicator[[i]]==sens_iterator
      
      if(sum(ind)==0){next}
      
      #If not mvrnorm or dirichlet, just update relevant value and leave remaining as they are 
      #Adjustment to mvrnorm or dirichlet is needed only if sum of parameters that need adjustment is < than length of elements in parameter
      if(!distributions[[i]] %in% c("mvrnorm","rdirichlet") | sum(ind)==length(temp_data[[i]])){
        output[[i]][ind] <- sens[[i]][ind] 
        
      } else if(distributions[[i]]=="mvrnorm"){
        #If distribution is mvrnorm, recalculate remaining values using conditional
        output[[i]] <- cond_mvn(temp_data[[i]], covariances[[i]], which(ind), sens[[i]][ind], full_output = TRUE)[[1]]
        
      } else { 
        #If distribution is dirichlet, recalculate remaining values using conditional
        output[[i]] <- cond_dirichlet(temp_data[[i]], which(ind), sens[[i]][ind], full_output = TRUE)
      }
    }
    
  } else{ #If using old approach of single parameters being changed simultaneously
    
    if ((any(!indicator %in% c(0,1)) & indicator_sens_binary == TRUE) | (!psa_ind %in% c(0,1)) | (!sens_ind %in% c(0,1)) ) {
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
    
    #if the parameter is out of the dsa/scenario specific iteration, or not in DSA/scenario, use PSA/normal. 
    for (it in 1:len_ind) {
      output[[it]] <-  if (indicator[[it]][1]==0 | sens_ind==F ) { 
        if (psa_ind==T & indicator_psa[it]==1) {
          psa[[it]]
        } else {
          base[[it]]
        }
      } else {#If active, use DSA if DSA and scenario if scenario
        sens[[it]]
      } 
    }
    
  }
  
  if (!is.null(names_out)) {
    names(output) <- names_out
  }
  
  return(as.list(output))
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
#' library(magrittr)
#'
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
#' @return No return value, adds event to `cur_evtlist` and integrates it with the main list for storage
#'
#' @importFrom stats setNames
#' 
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `new_event` with a list of one element, it's better to group them into a single `new_event` with a list of two elements.
#'
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {new_event(list("ae"=5))})


new_event <- function(evt){
  new_evt_name <- names(evt)
  new_evt <- setNames(unlist(evt),new_evt_name)
  if (!is.numeric(new_evt)) {
    stop("New event times are not all numeric, please review")
  }
  
  input_list_arm <- parent.frame()

  evtlist_temp <- list(cur_evtlist = c(input_list_arm$cur_evtlist,
                                  new_evt))

  if(input_list_arm$debug){ #only works correctly with create_if_null==TRUE, to be modified in later versions
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$sim,
                  "; Patient: ", input_list_arm$i,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(input_list_arm$curtime,3)
    )
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value,setNames(rep(Inf, length(new_evt_name)),new_evt_name))
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,new_evt)
      
    }else{
      dump_info <- list(
        list(prev_value = setNames(rep(Inf, length(new_evt_name)),new_evt_name),
             cur_value = new_evt
        )
      )
      names(dump_info) <- loc
      
      input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
      
    }
  }
  
  
  input_list_arm[["cur_evtlist"]] <- evtlist_temp$cur_evtlist

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
#' @return No return value, modifies/adds event to `cur_evtlist` and integrates it with the main list for storage
#'
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_event` with a list of one element, it's better to group them into a single `modify_event` with a list of two elements.
#'
#' This function does not evaluate sequentially.
#'
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {modify_event(list("os"=5))})

modify_event <- function(evt,create_if_null=TRUE){
  input_list_arm <- parent.frame()
  
  evt_unlist <- unlist(evt)
  if(!is.numeric(evt_unlist)){
      stop("Modify event times are not all numeric, please review")
  }
  names_evt <- names(evt)
  
  if(input_list_arm$debug){ #only works correctly with create_if_null==TRUE, to be modified in later versions
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$sim,
                  "; Patient: ", input_list_arm$i,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(input_list_arm$curtime,3)
                  )
    temp_cur <- input_list_arm$cur_evtlist[names_evt]
    isna_evt <- is.na(temp_cur)
    temp_cur[isna_evt] <- Inf
    names(temp_cur)[isna_evt] <- names_evt[isna_evt]
    
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value,temp_cur)
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,evt_unlist)
        
    }else{
    dump_info <- list(
      list(prev_value = temp_cur,
           cur_value = evt_unlist
           )
      )
    names(dump_info) <- loc
    
    input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
    
    }
  }
  
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
  

}





# Modify item in input list -------------------------------------------------------------------------------------------------------------------------------

#' Modify the value of existing items
#'
#' @param list_item A list of items and their values or expressions
#' 
#' @return No return value, modifies/adds item to the environment and integrates it with the main list for storage
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
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {modify_item(list("cost.it"=5))})


modify_item <- function(list_item){

  input_list_arm <- parent.frame()
  
  if(input_list_arm$debug){ 
    
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$sim,
                  "; Patient: ", input_list_arm$i,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(input_list_arm$curtime,3)
    )
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value, input_list_arm[names(list_item)])
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,list_item)
      
    }else{
    dump_info <- list(
      list(
        prev_value = input_list_arm[names(list_item)],
        cur_value = list_item
      )
    )
    names(dump_info) <- loc
    
    input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
    }
  }
  
  
  list2env(lapply(list_item, unname), parent.frame())
  
  if(input_list_arm$accum_backwards){
    l_temp <- as.list(rep(1,length(list_item)))
    names(l_temp) <- paste0(names(list_item),"_lastupdate",recycle0=TRUE)
    list2env(l_temp, parent.frame())
  }
  
  
}

# Modify item in input list  evaluating sequentially -------------------------------------------------------------------------------------------------------------------------------

#' Modify the value of existing items
#'
#' @param ... A list of items and their values or expressions. Will be evaluated sequentially (so one could have list(a= 1, b = a +2 ))
#'
#' @return No return value, modifies/adds items sequentially and deploys to the environment and with the main list for storage
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
#'  The function is different from modify_item in that this function evaluates sequentially the arguments within the list passed.
#'   This implies a slower performance relative to modify_item, but it can be more cleaner and convenient in certain instances.
#'
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {
#'   modify_item_seq(list(cost.idfs = 500, cost.tx = cost.idfs + 4000))
#'   })

modify_item_seq <- function(...){
  input_list_arm <- parent.frame()
  input_list <- as.list(substitute(...))[-1]
  list_out <- list()
  
  if(input_list_arm$debug){ 
    temp_dump <- mget(names(input_list),input_list_arm, ifnotfound = Inf)
  }
  
  for (inp in 1:length(input_list)) {
    temp_obj <- as.list(eval(input_list[[inp]], input_list_arm))
    names(temp_obj) <- names(input_list[inp])
    list2env(temp_obj, input_list_arm)
  }
  list_out <- mget(names(input_list),input_list_arm)
  
  if(input_list_arm$debug){ 
    
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$sim,
                  "; Patient: ", input_list_arm$i,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(input_list_arm$curtime,3)
    )
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value, temp_dump)
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,list_out)
      
    }else{
    dump_info <- list(
      list(prev_value = temp_dump,
           cur_value = list_out
      )
    )
    names(dump_info) <-loc
    input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
    }
  }
  
  
  
  
  # list2env(list_out,envir = parent.frame())

  if(input_list_arm$accum_backwards){
    l_temp <- as.list(rep(1,length(input_list)))
    names(l_temp) <- paste0(names(input_list),"_lastupdate",recycle0=TRUE)
    list2env(l_temp, parent.frame())
  }
  

}

# Add_reactevt -------------------------------------------------------------------------------------------------------------------------------------------

#' Define the modifications to other events, costs, utilities, or other items affected by the occurrence of the event
#'
#' @param .data Existing data for event reactions
#' @param name_evt Name of the event for which reactions are defined.
#' @param input Expressions that define what happens at the event, using functions as defined in the Details section
#'
#' @return A named list with the event name, and inside it the substituted expression saved for later evaluation
#'
#' @export
#'
#' @details
#' There are a series of objects that can be used in this context to help define the event reactions.
#'
#' The following functions may be used to define event reactions within this `add_reactevt()` function:
#' `modify_item()` | Adds & Modifies items/flags/variables for future events (does not consider sequential)
#' `modify_item_seq()` | Adds & Modifies items/flags/variables for future events in a sequential manner
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
#' The user can use `extract_from_reactions` function on the output to obtain a data.frame with all the relationships defined in the reactions in the model.
#'
#' @examples
#' add_reactevt(name_evt = "start",input = {})
#' add_reactevt(name_evt = "idfs",input = {modify_item(list("fl.idfs"= 0))})

add_reactevt <- function(.data=NULL,name_evt,input){

  if (length(name_evt)>1 | !is.character(name_evt) | any(nchar(name_evt)<2)) {
    stop("name_evt argument in add_reactevt should be a single string with at least 2 characters")
  }
  
  data_list <- .data
  
  
  input_sub <- substitute(input)
  ret_sub <- substitute(return(environment()))
  new_sub <- as.call(c(substitute(`{`), as.list(input_sub)[-1], ret_sub))

  evt_r <- list(list(react=new_sub))
  
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
#' add_tte(arm="int",evts = c("start","ttot","idfs","os"),
#' input={
#' start <- 0
#' idfs <- draw_tte(1,'lnorm',coef1=2, coef2=0.5)
#' ttot <- min(draw_tte(1,'lnorm',coef1=1, coef2=4),idfs)
#' os <- draw_tte(1,'lnorm',coef1=0.8, coef2=0.2)
#' })
#'
add_tte <- function(.data=NULL,arm, evts, other_inp = NULL,input){
  data_list <- .data

  for (arm_i in arm) {

    if (!is.character(other_inp) & !is.null(other_inp)) {
      stop("other_inp argument is required to be a character vector or be set to NULL")
    }

    
    if (length(evts)<2 | !is.character(evts) | any(nchar(evts)<2)) {
      stop("evts argument in add_tte should be a string vector with at least 2 characters each. At least two events should be defined (e.g., start and end).")
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

# Perform luck adjustment when using conditional quantile -------------------------------------------------------------------------------------------------------------------------------

#' Perform luck adjustment
#'
#' @param prevsurv Value of the previous survival
#' @param cursurv Value of the current survival
#' @param luck Luck used to be adjusted (number between 0 and 1)
#' @param condq Conditional quantile approach or standard approach
#'
#' @export
#' 
#' @return Adjusted luck number between 0 and 1
#'
#' @details
#' This function performs the luck adjustment automatically for the user, returning the adjusted luck number.
#' Luck is interpreted in the same fashion as is standard in R (higher luck, higher time to event).
#' 
#' Note that if TTE is predicted using a conditional quantile function (e.g., conditional gompertz, conditional quantile weibull...) `prevsurv` and `cursurv`
#' are the unconditional survival using the "previous" parametrization but at the previous time for `presurv` and at the current time for `cursurv`.
#' For other distributions, `presurv` is the survival up to current time using the previous parametrization, and `cursurv` 
#' is the survival up to current time using the current parametrization.
#' 
#' Note that the advantage of the conditional quantile function is that it does not need the new parametrization to update the luck, 
#' which makes this approach computationally more efficient.
#' This function can also work with vectors, which could allow to update multiple lucks in a single approach, and it can preserve names
#'
#' @examples
#' luck_adj(prevsurv = 0.8,
#'  cursurv = 0.7,
#'  luck = 0.5,
#'  condq = TRUE)
#'  
#' luck_adj(prevsurv = c(1,0.8,0.7),
#'  cursurv = c(0.7,0.6,0.5),
#'  luck = setNames(c(0.5,0.6,0.7),c("A","B","C")),
#'  condq = TRUE)
#'  
#' luck_adj(prevsurv = 0.8,
#'  cursurv = 0.7,
#'  luck = 0.5,
#'  condq = FALSE) #different results
#' 
#' #Unconditional approach, timepoint of change is 25,
#' # parameter goes from 0.02 at time 10 to 0.025 to 0.015 at time 25,
#' #  starting luck is 0.37
#' new_luck <- luck_adj(prevsurv = 1 - pweibull(q=10,3,1/0.02),
#'  cursurv = 1 - pweibull(q=10,3,1/0.025),
#'  luck = 0.37,
#'  condq = FALSE) #time 10 change
#'  
#' new_luck <- luck_adj(prevsurv = 1 - pweibull(q=25,3,1/0.025),
#'  cursurv = 1 - pweibull(q=25,3,1/0.015),
#'  luck = new_luck,
#'  condq = FALSE) #time 25 change
#'  
#' qweibull(new_luck, 3, 1/0.015) #final TTE 
#' 
#' #Conditional quantile approach 
#' new_luck <- luck_adj(prevsurv = 1-pweibull(q=0,3,1/0.02),
#'                       cursurv = 1- pweibull(q=10,3,1/0.02),
#'                       luck = 0.37,
#'                       condq = TRUE) #time 10 change, previous time is 0 so prevsurv will be 1
#' 
#' new_luck <- luck_adj(prevsurv = 1-pweibull(q=10,3,1/0.025),
#'                       cursurv = 1- pweibull(q=25,3,1/0.025),
#'                       luck = new_luck,
#'                       condq = TRUE) #time 25 change
#' 
#' qcond_weibull(rnd = new_luck,
#'                      shape = 3,
#'                      scale = 1/0.015,
#'                      lower_bound = 25) + 25 #final TTE

luck_adj <- function(prevsurv,cursurv,luck,condq=TRUE){
  
  #If length is 1, go for faster approach, otherwise use vectorial approach
  if(length(prevsurv)==1){
    if(condq==TRUE){
      adj_luck <-  max(1e-5, min(0.99999, if(prevsurv != 0){ 1 - ((1 - luck) * (prevsurv/cursurv)) } else{ luck }))
    }else{
      adj_luck <-  max(1e-5, min(0.99999, if(prevsurv != 0){ 1 - ((1 - luck) * (cursurv/prevsurv)) } else{ luck }))
    }
  } else{
    if(condq==TRUE){
      adj_luck <-  pmax(1e-5, pmin(0.99999, ifelse(prevsurv != 0, 1 - ((1 - luck) * (prevsurv/cursurv)), luck )))
    }else{
      adj_luck <-   pmax(1e-5, pmin(0.99999, ifelse(prevsurv != 0, 1 - ((1 - luck) * (cursurv/prevsurv)), luck )))
    }
  }
  
  #preserve names if existing
  names_luck <- names(luck)
  if(!is.null(names_luck)){
    names(adj_luck) <-names_luck
  }
  
  return(adj_luck)
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
#' @examples 
#' disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' 
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
#' @examples 
#' disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' 
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
#' @examples
#' disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500)
#' 
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
#' @examples 
#' disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500)
#' 
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
#' @examples 
#' disc_cycle(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' 
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
#' @examples 
#' disc_cycle_v(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' 
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

# Model Events Interactions -------------------------------------------------------

#' Extract all items and events and their interactions from the event reactions list
#'
#' @param reactions list generated through `add_reactevt`
#'
#' @return A data.frame with the relevant item/event, the event where it's assigned,
#'  and whether it's contained within a conditional statement
#'
#' @examples 
#' evt_react_list2 <-
#'   add_reactevt(name_evt = "sick",
#'                input = {modify_item(list(a=1+5/3))
#'                  assign("W", 5 + 3 / 6 )
#'                  x[5] <- 18
#'                  for(i in 1:5){
#'                    assign(paste0("x_",i),5+3)
#'                  }
#'                  if(j == TRUE){
#'                    y[["w"]] <- 612-31+3
#'                  }#'                
#'                  q_default <- 0
#'                  c_default <- 0
#'                  curtime   <- Inf
#'                })
#'    
#'  extract_from_reactions(evt_react_list2)
#' 
#'
#' @export
#' 
extract_from_reactions <- function(reactions){
  data.table::rbindlist(lapply(purrr::map(reactions,"react"), function(x){
    out <- ast_as_list(x)
    extract_elements_from_list(out)
  }),
  idcol="event")
}


#' Transform a substituted expression to its Abstract Syntax Tree (AST) as a list
#'
#' @param ee Substituted expression
#'
#' @return Nested list with the Abstract Syntax Tree (AST)
#'
#' @export
#'
#' @examples
#' expr <- substitute({
#' 
#' a <- sum(5+7)
#' 
#' modify_item(list(afsa=ifelse(TRUE,"asda",NULL)))
#' 
#' modify_item_seq(list(
#'   
#'   o_other_q_gold1 = if(gold == 1) { utility } else { 0 },
#'   
#'   o_other_q_gold2 = if(gold == 2) { utility } else { 0 },
#'   
#'   o_other_q_gold3 = if(gold == 3) { utility } else { 0 },
#'   
#'   o_other_q_gold4 = if(gold == 4) { utility } else { 0 },
#'   
#'   o_other_q_on_dup = if(on_dup) { utility } else { 0 }
#'  
#' ))
#' 
#' if(a==1){
#'   modify_item(list(a=list(6+b)))
#'   
#'   modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#' } else{
#'   modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#'   if(a>6){
#'     modify_item(list(a=8))
#'   }
#'   
#' }
#' 
#' 
#' if (sel_resp_incl == 1 & on_dup == 1) {
#'   
#'   modify_event(list(e_response = curtime, z = 6))
#'   
#' }
#' 
#' })
#' 
#' 
#' out <- ast_as_list(expr)
#' 
ast_as_list <- function(ee) {
  purrr::map_if(as.list(ee), is.call, ast_as_list)
  }


#' Extracts items and events by looking into assignments, modify_item, modify_item_seq, modify_event and new_event
#'
#' @param node Relevant node within the nested AST list
#' @param conditional_flag Boolean whether the statement is contained within a conditional statement
#'
#' @return A data.frame with the relevant item/event, the event where it's assigned,
#'  and whether it's contained within a conditional statement
#'  
#' @examples
#' expr <- substitute({
#' 
#' a <- sum(5+7)
#' 
#' modify_item(list(afsa=ifelse(TRUE,"asda",NULL)))
#' 
#' modify_item_seq(list(
#'   
#'   o_other_q_gold1 = if(gold == 1) { utility } else { 0 },
#'   
#'   o_other_q_gold2 = if(gold == 2) { utility } else { 0 },
#'   
#'   o_other_q_gold3 = if(gold == 3) { utility } else { 0 },
#'   
#'   o_other_q_gold4 = if(gold == 4) { utility } else { 0 },
#'   
#'   o_other_q_on_dup = if(on_dup) { utility } else { 0 }
#'  
#' ))
#' 
#' if(a==1){
#'   modify_item(list(a=list(6+b)))
#'   
#'   modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#' } else{
#'   modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#'   if(a>6){
#'     modify_item(list(a=8))
#'   }
#'   
#' }
#' 
#' 
#' if (sel_resp_incl == 1 & on_dup == 1) {
#'   
#'   modify_event(list(e_response = curtime, z = 6))
#'   
#' }
#' 
#' })
#' 
#' 
#' out <- ast_as_list(expr)
#' 
#' results <- extract_elements_from_list(out)
#' 
#' 
#' @export
#' 
extract_elements_from_list <- function(node, conditional_flag = FALSE) {
  results <- data.frame(
    name = character(),
    type = character(),
    conditional_flag = logical(),
    definition = character(),
    stringsAsFactors = FALSE
  )
  
  if (is.list(node)) {
    func_name <- if (!is.null(node[[1]])) as.character(node[[1]]) else NULL
    
    # Case 1: modify_* / new_event
    if (!is.null(func_name) && any(func_name %in% c("modify_item_seq", "modify_item", "modify_event", "new_event"))) {
      type <- if (any(func_name %in% c("modify_item_seq", "modify_item"))) "item" else "event"
      list_expr <- node[[2]]
      if (is.list(list_expr)) {
        definition <- unlist(extract_defs(node), recursive = TRUE)
        definition <- definition[!names(definition) == ""]
        results_temp <- data.frame(
          name = names(definition),
          type = type,
          conditional_flag = conditional_flag,
          definition = definition,
          stringsAsFactors = FALSE
        )
        results <- rbind(results, results_temp)
      }
    }
    
    # Case 2: assignment via <- or =
    if (!is.null(func_name) && func_name %in% c("<-", "=")) {

      lhs <- node[[2]]
      rhs <- node[[3]]
      # Detect and extract from complex LHS expressions
      if (is.symbol(lhs)) {
        target_name <- as.character(lhs)
      } else if (is.list(lhs)) {
        target_name <- clean_output(expr_from_list(lhs))
      } else {
        target_name <- NA
      }

      if (!is.na(target_name)) {
        results <- rbind(results, data.frame(
          name = target_name,
          type = "item",
          conditional_flag = conditional_flag,
          definition = clean_output(expr_from_list(rhs)),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Case 3: assignment via assign()
    if (!is.null(func_name) && func_name == "assign") {
      if (length(node) >= 3) {
        varname <- node[[2]]
        if (is.character(varname)) {
          varname <- as.character(varname)
        } else if (is.list(varname)) {
          varname <- clean_output(expr_from_list(varname))
        }
        
        value_expr <- node[[3]]
        if (is.character(varname)) {
          results <- rbind(results, data.frame(
            name = varname,
            type = "item",
            conditional_flag = conditional_flag,
            definition = clean_output(expr_from_list(value_expr)),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Case 4: if statement → propagate conditional flag
    if (!is.null(func_name) && func_name == "if") {
      conditional_flag <- TRUE
    }
  }
  
  # Recursively walk children
  if (is.list(node)) {
    for (child in node) {
      results <- rbind(results, extract_elements_from_list(child, conditional_flag))
    }
  }
  
  results <- results[!(is.na(results$name) | results$name == "" | is.na(results$definition)), ]
  results
}


#' Loop to extract an expression from a list
#'
#' @param lst sublist from the AST as list
#'
#' @return Reconstructed expression as a character
#'
#' @noRd
#' 
expr_from_list <- function(lst) {
  if(is.null(lst)){
    return(lst)
  } else if (is.atomic(lst)) {
    return(lst)
  } else if (is.symbol(lst)) {
    return(as.name(lst))
  } else if (is.list(lst) && length(lst) == 1) {
    return(expr_from_list(lst[[1]]))
  } else {
    func_name <- deparse(expr_from_list(lst[[1]]))
    args <- lapply(lst[-1], expr_from_list)
    
    # Handle special cases (e.g., operators)
    return(do.call(call, c(func_name, args),quote = TRUE))
  }
}


#' Clean output from generated expression
#'
#' @param called Generated expression from `expr_from_list`
#'
#' @return Cleaned character expression 
#'
#' @noRd
#' 
clean_output <- function(called){
  gsub("\`","",
       gsub("    ","",
            gsub("\"","'",
            paste0(
              deparse(called
              ),collapse="")
            )
      )
  )
}

#' Extract relevant definitions from sublist
#'
#' @param lst sublist of the AST where to obtain definitions from
#'
#' @return Vector of cleaned character expressions definitions found in the sublist 
#'
#' @noRd
#' 
extract_defs <- function(lst){
  if(purrr::pluck_depth(lst)<=1 & is.null(names(lst))){
    NULL
  }else{
    if(length(lst[!is.null(names(lst))&names(lst)!=""])==0){
      lapply(lst, extract_defs)
    } else{
      lapply(lst[!is.null(names(lst)) &names(lst)!=""], function(x){
        if(is.null(x)){"NULL"}else{
          clean_output(expr_from_list(x))              
          
        }
      }
      )
    }
  }
}



#' Extract assignments from expression for debug mode
#'
#' @param expr Expression of type language 
#'
#' @return Vector of cleaned character objects that are assigned. It will ignore objects created dynamically (e.g., by pasting in a loop)
#'
#' @noRd
#'  
#' @examples
#' 
#'  expr <- substitute({
#'    q_default <- if (fl.idfs==1) {
#'      util.idfs.ontx * fl.idfs.ontx + (1-fl.idfs.ontx) * (1-fl.idfs.ontx) 
#'    } else if (fl.idfs==0 & fl.mbcs==0) {
#'      util.remission * fl.remission + fl.recurrence*util.recurrence
#'    } else if (fl.mbcs==1) {
#'      util.mbc.progression.mbc * fl.mbcs.progression.mbc + (1-fl.mbcs.progression.mbc)*util.mbc.pps
#'    }
#'    c_default <- cost.recurrence * fl.recurrence
#'    fl.recurrence <- 1
#'    fl.remission <- 0
#'    fl.mbcs <- 1
#'    fl.mbcs.progression.mbc <- 1 #ad-hoc for plot
#'    if(a){c<-5}else{d<-6}
#'    if(e){assign("t",6)}else{j<-6}
#'    for(i in 1:10){
#'      assign(paste0("a_",i),1)
#'    }
#'    a = 5
#'    d <- b <- c <- 5
#'  })
#'  
#'  extract_assignment_targets(expr)
#' 
#' 
extract_assignment_targets <- function(expr) {
  assigned <- character()
  
  walk_node <- function(node) {
    if (is.call(node)) {
      fname <- as.character(node[[1]])
      
      # Handle <- and = assignment
      if (fname %in% c("<-", "=")) {
        lhs <- node[[2]]
        if (is.symbol(lhs)) {
          assigned <<- c(assigned, as.character(lhs))
        } else if (is.call(lhs) && lhs[[1]] == as.name("$")) {
          # x$y <- ...  => get "x"
          assigned <<- c(assigned, as.character(lhs[[2]]))
        } else if (is.call(lhs) && lhs[[1]] == as.name("[[")) {
          # x[["y"]] <- ... => get "x"
          assigned <<- c(assigned, as.character(lhs[[2]]))
        } else if (is.call(lhs) && lhs[[1]] == as.name("[")) {
          assigned <<- c(assigned, as.character(lhs[[2]]))
        }
        walk_node(node[[3]])  # Recurse on RHS
      }
      
      # Handle assign("var", value)
      else if (fname == "assign" && length(node) >= 2) {
        varname <- node[[2]]
        if (is.character(varname)) {
          assigned <<- c(assigned, varname)
        } else if (is.symbol(varname)) {
          assigned <<- c(assigned, as.character(varname))
        }
        # Recurse on value
        if (length(node) >= 3) walk_node(node[[3]])
      }
      
      # Recurse into all sub-calls
      for (arg in as.list(node)[-1]) {
        walk_node(arg)
      }
    } else if (is.expression(node) || is.call(node)) {
      for (item in as.list(node)) walk_node(item)
    }
  }
  
  walk_node(expr)
  unique(assigned)
}