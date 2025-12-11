
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
#' @param deploy_env Boolean, if TRUE will deploy all objects in the environment where the function is called for. Must be active if using add_item (and FALSE if a list must be returned)
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
#'              distributions = list("rnorm","rnorm"),
#'              deploy_env = FALSE
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
#'              covariances = list(0.1,0.1,matrix(c(1,0.1,0.1,1),2,2)),
#'              deploy_env = FALSE
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
                         covariances = NULL,
                         deploy_env = TRUE
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
  
  if(deploy_env){
    list2env(as.list(output), parent.frame())
  } else{
    return(as.list(output))
  }
} 


# Add item/parameter to list --------------------------------------------------------

#' Define or append model inputs
#'
#' Build a single `{}` expression that defines inputs for a simulation.
#' - Named args in `...` become assignments (`name <- expr`), e.g., `add_item(a=5)`
#' - Unnamed args are inserted raw/unevaluated. If an unnamed arg is a `{}` block,
#'   its statements are spliced (flattened). `add_item(pick_val_v(...))`
#' - Works with magrittr pipes: a leading `.` (the LHS) is resolved to its value;
#'   if that value is a `{}` block (or list of expressions), it becomes the
#'   starting block.
#' - input argument can be used to handle alternative `add_item2` method,
#'   e.g. `add_item(input = {a <- 5})`
#'
#' @param ... Unevaluated arguments. Named → `name <- expr`; unnamed → raw expr.
#' @param .data Optional named argument: an existing `{}` block (or list of
#'   expressions) to start from.
#' @param input Optional unevaluated expression or `{}` block to splice in.
#'
#' @return A single `{}` call (language object) ready for `load_inputs()`.
#' 
#' @examples
#' library(magrittr)
#'
#' add_item(input = {fl.idfs <-  0})
#' add_item(input = {
#'  util_idfs <- if(psa_bool){rnorm(1,0.8,0.2)} else{0.8}
#'  util.mbc <- 0.6
#'  cost_idfs <- 2500})
#' common_inputs <- add_item(input = {
#' pick_val_v(
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
#'   names_out = l_statics[["parameter_name"]],
#'   deploy_env = TRUE #Note this option must be active if using it at add_item2
#' )
#' }
#' ) 
#' 
#' @export
add_item <- function(..., .data = NULL, input) {
  mc   <- match.call(expand.dots = FALSE)
  dots <- mc$...
  
  # Helper: coerce various starters into a { } block list
  as_block_list <- function(x) {
    if (is.null(x)) {
      list(as.name("{"))
    } else if (is.call(x) && identical(x[[1L]], as.name("{"))) {
      as.list(x)
    } else if (is.list(x)) {
      c(list(as.name("{")), x)
    } else {
      c(list(as.name("{")), list(x))
    }
  }
  
  # 1) Start from .data if provided
  block_elems <- as_block_list(.data)
  built <- list()
  
  splice_or_keep <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("{"))) {
      as.list(expr)[-1L]
    } else {
      list(expr)
    }
  }
  
  # 2) If .data is NULL and first unnamed ... is DOT (.) from magrittr, resolve it
  if (is.null(.data) && !is.null(dots) && length(dots) > 0) {
    dn <- names(dots)
    first_is_named <- !is.null(dn) && nzchar(dn[1L])
    if (!first_is_named) {
      first_expr <- dots[[1L]]
      if (is.symbol(first_expr) && identical(first_expr, as.name("."))) {
        # Evaluate DOT to get LHS value from the calling frame
        lhs <- eval(first_expr, parent.frame())
        block_elems <- as_block_list(lhs)
        dots <- if (length(dots) > 1L) dots[-1L] else NULL
      } else if (is.call(first_expr) && identical(first_expr[[1L]], as.name("{"))) {
        # Promote a literal { } block as the starter
        block_elems <- as.list(first_expr)
        dots <- if (length(dots) > 1L) dots[-1L] else NULL
      }
    }
  }
  
  # 3) Process remaining ... : named -> assignment; unnamed -> raw/splice
  if (!is.null(dots)) {
    dn <- names(dots)
    for (i in seq_along(dots)) {
      nm   <- if (length(dn)) dn[[i]] else NULL
      expr <- dots[[i]]
      if (!is.null(nm) && nzchar(nm)) {
        built[[length(built) + 1L]] <- call("<-", as.name(nm), expr)
      } else {
        built <- c(built, splice_or_keep(expr))
      }
    }
  }
  
  # 4) Optional input=
  if (!missing(input)) {
    input_sub <- substitute(input)
    built <- c(built, splice_or_keep(input_sub))
  }
  
  # 5) Return a proper { ... } call
  as.call(c(block_elems, built))
}

#' Define parameters that may be used in model calculations (uses expressions)
#'
#' @param .data Existing data
#' @param input Items to define for the simulation as an expression (i.e., using {})
#'
#' @return A substituted expression to be evaluated by engine 
#'
#' @importFrom lifecycle deprecate_stop
#'
#' @export
#'
#' @details
#' DEPRECATED (old description): The functions to add/modify events/inputs use named vectors or lists. If chaining together add_item2, it will just append the expressions together in the order established.
#'
#' If using `pick_val_v`, note it should be used with the `deploy_env = TRUE` argument so that add_item2 process it correctly.
#'
add_item2 <- function(.data=NULL,input){
  lifecycle::deprecate_stop("2.0.0", "add_item2()", "add_item()")
  
  data_list <- .data
  
  list_item <- substitute(input)
  
  if (is.null(data_list)) {
    data_list <- list_item
  } else{
    data_list <- as.call(c(substitute(`{`), as.list(data_list)[-1], as.list(list_item)[-1]))
  }
  
  return(data_list)
}

# Modify item in input list -------------------------------------------------------------------------------------------------------------------------------

#' Modify the value of existing items
#'
#' @param list_item A list of items and their values or expressions
#' 
#' @return No return value, modifies/adds item to the environment and integrates it with the main list for storage
#'
#' @importFrom lifecycle deprecate_stop
#'
#' @export
#'
#' @details
#' DEPRECATED (old description): The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_item` with a list of one element, it's better to group them into a single `modify_item` with a list of two elements.
#' 
#' Note that `modify_item` nor `modify_item_seq` can work on subelements (e.g.,
#' `modify_item(list(obj$item = 5))` will not work as intended, for that is better
#' to assign directly using the expression approach, so `obj$item <- 5`).
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
  lifecycle::deprecate_stop("2.0.0", "modify_item()",
                            details = "This function is deprecated, as it's no longer needed to run the simulation,
                            not even when backwards = TRUE. Just write a standard code (see examples).")
  input_list_arm <- parent.frame()

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
#' @importFrom lifecycle deprecate_stop
#'
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_item` with a list of one element, it's better to group them into a single `modify_item` with a list of two elements.
#'
#' Note that `modify_item` nor `modify_item_seq` can work on subelements (e.g.,
#' `modify_item_seq(list(obj$item = 5))` will not work as intended, for that is better
#' to assign directly using the expression approach, so `obj$item <- 5`).
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
  lifecycle::deprecate_stop("2.0.0", "modify_item_seq()",
                            details = "This function is deprecated, as it's no longer needed to run the simulation,
                            not even when backwards = TRUE. Just write a standard code (see examples).")
  
  input_list_arm <- parent.frame()
  input_list <- as.list(substitute(...))[-1]
  list_out <- list()

  for (inp in 1:length(input_list)) {
    temp_obj <- as.list(eval(input_list[[inp]], input_list_arm))
    names(temp_obj) <- names(input_list[inp])
    list2env(temp_obj, input_list_arm)
  }
  list_out <- mget(names(input_list),input_list_arm)
  
  # list2env(list_out,envir = parent.frame())

  if(input_list_arm$accum_backwards){
    l_temp <- as.list(rep(1,length(input_list)))
    names(l_temp) <- paste0(names(input_list),"_lastupdate",recycle0=TRUE)
    list2env(l_temp, parent.frame())
  }
  

}


# event queueing ----------------------------------------------------------

#' Create a New Event Queue
#'
#' Initializes a new event queue with the specified priority order of event names.
#'
#' @param priority_order A character vector of event names sorted by decreasing importance.
#'
#' @return An external pointer to the new event queue.
#' @export
queue_create <- function(priority_order) {
  queue_create_cpp(priority_order)
}

#' Add events to the queue for a patient
#'
#' Adds one or more events for a given patient to the queue.
#'
#' @param events A named numeric vector. Names are event types, values are event times. It can also handle lists instead of named vectors (at a small computational cost).
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#'
#' @return NULL (invisible). Modifies the queue in-place.
#' @export
#' 
#' @importFrom stats setNames
#'
#' @details
#' The functions to add/modify events/inputs use named vectors or lists. Whenever several inputs/events are added or modified,
#'  it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `new_event` with a list of one element, it's better to group them into a single `new_event` with a list of two elements.
#' 
#' While multiple events can be added, they must be named differently. If the same event is added multiple times at once, only the last occurrence will be kept
#'  (only one event per event type in the queue of events yet to occur). If an event occurs, then a new one with the same name can be set.
#'
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {new_event(c("ae"=5))})
new_event <- function(events, ptr, patient_id) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  if (missing(patient_id)) patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
  if(is.list(events)){
    events <- unlist(events)
  }
  
  input_list_arm <- parent.frame()
  if(is.null(input_list_arm$debug)){input_list_arm$debug <- FALSE}
  if(input_list_arm$debug){
    new_evt_name <- names(events)
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$simulation,
                  "; Patient: ", patient_id,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(ifelse(is.null(input_list_arm$curtime),NA,input_list_arm$curtime),3)
    )
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value,setNames(rep(Inf, length(new_evt_name)),new_evt_name))
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,events)

    }else{
      dump_info <- list(
        list(prev_value = setNames(rep(Inf, length(new_evt_name)),new_evt_name),
             cur_value = events
        )
      )
      names(dump_info) <- loc

      input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)

    }
  }
  
  new_event_cpp(ptr, patient_id, events)
}

#' Get the next events in the queue
#'
#' Retrieves the next `n` events (without removing them).
#'
#' @param n Number of events to retrieve. Default is 1.
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#'
#' @return A list of events, each with `patient_id`, `event_name`, and `time`.
#' @export
next_event <- function(n = 1, ptr) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  next_event_cpp(ptr, n)
}

#' Get the next events in the queue for a specific patient
#'
#' Retrieves the next `n` events (without removing them).
#'
#' @param n Number of events to retrieve. Default is 1.
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#'
#' @return A list of events, each with `patient_id`, `event_name`, and `time`.
#' @export
next_event_pt <- function(n = 1, ptr, patient_id) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  if (missing(patient_id)) patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
  
  next_event_pt_cpp(ptr, patient_id, n)
}


#' Remove the next event from the queue
#'
#' Removes the next scheduled event from the queue. Not needed by user.
#'
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#'
#' @return NULL (invisible). Modifies the queue in-place.
pop_event <- function(ptr) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  pop_event_cpp(ptr)
}

#' Pop and return the next event
#'
#' Removes the next event from the queue and returns its details. Not needed by user.
#'
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#'
#' @return A named list with `patient_id`, `event_name`, and `time`.
pop_and_return_event <- function(ptr) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  pop_and_return_event_cpp(ptr)
}

#' Remove events for a patient
#'
#' Removes one or more events from the queue for the given patient.
#'
#' @param events A character vector of event names to remove. It can also handle lists instead of named vectors (at a small computational cost).
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#'
#' @return NULL (invisible). Modifies the queue in-place.
#' @export
remove_event <- function(events, ptr, patient_id) {
  if (missing(ptr)) {ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)}
  if (missing(patient_id)) {patient_id <- get("i", envir = parent.frame(), inherits = TRUE)}
  if(is.list(events)){
    events <- unlist(events)
  }
  
  input_list_arm <- parent.frame()
  if(is.null(input_list_arm$debug)){input_list_arm$debug <- FALSE}
  if(input_list_arm$debug){
    events_values <- setNames(get_event(events, ptr, patient_id),events)
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$simulation,
                  "; Patient: ", patient_id,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(ifelse(is.null(input_list_arm$curtime),NA,input_list_arm$curtime),3)
    )
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value,events_values)
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,setNames(rep(Inf, length(events)),events))
      
    }else{
      dump_info <- list(
        list(prev_value = events_values,
             cur_value = setNames(rep(Inf, length(events)),events)
        )
      )
      names(dump_info) <- loc
      
      input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
      
    }
  }
  
  remove_event_cpp(ptr, patient_id, events)
}

#' Modify or add events for a patient
#'
#' Modifies existing event times, or adds new events if `create_if_missing` is TRUE.
#'
#' @param events A named numeric vector with event names and new event times. It can also handle lists instead of named vectors (at a small computational cost).
#' @param create_if_missing Logical, whether to create events if they do not exist.
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#'
#' @return NULL (invisible). Modifies the queue in-place.
#' @export
#'
#' @details
#' The functions to add/modify events/inputs use named vectors or lists. Whenever several inputs/events are added or modified, it's recommended to group them within one function, as it reduces the computation cost.
#' So rather than use two `modify_event` with a list of one element, it's better to group them into a single `modify_event` with a list of two elements.
#'
#' This function does not evaluate sequentially.
#' 
#' While multiple events can be added, they must be named differently. If the same event is added multiple times at once, only the last occurrence will be kept
#'  (only one event per event type in the queue of events yet to occur). If an event occurs, then a new one with the same name can be set.
#'
#' This function is intended to be used only within the `add_reactevt` function in its `input` parameter and should not be run elsewhere or it will return an error.
#'
#' @examples
#' add_reactevt(name_evt = "idfs",input = {modify_event(c("os"=5))})
modify_event <- function(events, create_if_missing = TRUE, ptr, patient_id) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  if (missing(patient_id)) patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
  if(is.list(events)){
    events <- unlist(events)
  }
  
  input_list_arm <- parent.frame()
  if(is.null(input_list_arm$debug)){input_list_arm$debug <- FALSE}
  if(input_list_arm$debug){ 
    loc <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                  "; Sim: ", input_list_arm$simulation,
                  "; Patient: ", patient_id,
                  "; Arm: ", input_list_arm$arm,
                  "; Event: ", input_list_arm$evt,
                  "; Time: ", round(input_list_arm$curtime,3)
    )
    
    temp_cur <- c()
    for (idx in 1:length(events)) {
      temp_cur[idx] <- tryCatch(get_event(events[idx],ptr,patient_id), error = function(e) NA)
    }
    names(temp_cur) <- names(events)
    
    if(!is.null(input_list_arm$log_list[[loc]])){
      input_list_arm$log_list[[loc]]$prev_value <- c(input_list_arm$log_list[[loc]]$prev_value,temp_cur)
      input_list_arm$log_list[[loc]]$cur_value <- c(input_list_arm$log_list[[loc]]$cur_value,events)
      
    }else{
      dump_info <- list(
        list(prev_value = temp_cur,
             cur_value = events
        )
      )
      names(dump_info) <- loc
      
      input_list_arm$log_list <- c(input_list_arm$log_list, dump_info)
      
    }
  }
  
  modify_event_cpp(ptr, patient_id, events, create_if_missing)
}

#' Check if the event queue is empty
#'
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param exclude_inf Logical, whether to exclude events with Inf time. Default is FALSE.
#'
#' @return Logical, TRUE if the queue is empty, FALSE otherwise.
#' @export
queue_empty <- function(ptr, exclude_inf = FALSE) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  queue_empty_cpp(ptr, exclude_inf)
}

#' Get the Size of the Event Queue
#'
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param exclude_inf Logical, whether to exclude events with Inf time. Default is FALSE.
#'
#' @return An integer indicating the number of events in the queue.
#' @export
queue_size <- function(ptr, exclude_inf = FALSE) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  queue_size_cpp(ptr, exclude_inf)
}

#' Check if a patient has a specific event in the queue
#'
#' @param event_name Character string, the name of the event.
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#' @param exclude_inf Logical, whether to exclude events with Inf time. Default is FALSE.
#'
#' @return Logical, TRUE if the event exists for the patient (optionally excluding Inf), FALSE otherwise.
#' @export
has_event <- function(event_name, ptr, patient_id, exclude_inf = FALSE) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  if (missing(patient_id)) patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
  has_event_cpp(ptr, patient_id, event_name, exclude_inf)
}


#' Get a specific event time
#'
#' @param event_name Character string, the name of the event.
#' @param ptr The event queue pointer. Defaults to `cur_evtlist`.
#' @param patient_id The patient ID. Defaults to `i`.
#'
#' @return Numeric, time of event for patient
#' @export
get_event <- function(event_name, ptr , patient_id ) {
  if (missing(ptr)) ptr <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  if (missing(patient_id)) patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
  get_event_cpp(ptr, patient_id, event_name)
}


# discrete resource -------------------------------------------------------

#' Create a discrete resource
#' 
#' Creates a discrete resource management system for discrete event simulations.
#' This system manages a fixed number of identical resource units that can be
#' blocked (used) by patients and maintains a priority queue for waiting patients.
#' 
#' @param n Integer. The total capacity of the resource (must be >= 1).
#' 
#' @return An environment with methods for resource management.
#' 
#' @details
#' The returned environment has the following methods:
#' \itemize{
#'   \item \code{size()}: Returns the total capacity
#'   \item \code{queue_size()}: Returns the number of patients in queue
#'   \item \code{n_free()}: Returns the number of free resource units
#'   \item \code{patients_using()}: Vector of patient IDs currently using the resource
#'   \item \code{patients_using_times()}: Vector of start times for patients using the resource
#'   \item \code{queue_start_times()}: Vector of queue start times parallel to queue order
#'   \item \code{queue_priorities()}: Vector of priorities parallel to queue order
#'   \item \code{queue_info(n)}: Data.frame with patient_id, priority, start_time for queue
#'   \item \code{is_patient_in_queue(patient_id)}: Check if patient is in queue
#'   \item \code{is_patient_using(patient_id)}: Check if patient is using resource
#'   \item \code{attempt_block(patient_id, priority, start_time)}: Attempt to block a resource unit
#'   \item \code{attempt_free(patient_id, remove_all)}: Free a resource unit
#'   \item \code{attempt_free_if_using(patient_id, remove_all)}: Free only if patient is using
#'   \item \code{next_patient_in_line(n)}: Get next n patients in queue
#'   \item \code{modify_priority(patient_id, new_priority)}: Modify patient priority in queue
#'   \item \code{add_resource(n)}: Add n resource units to total capacity
#'   \item \code{remove_resource(n, current_time)}: Remove n resource units from total capacity
#' }
#' 
#' @examples
#' # Create a resource with 3 units
#' beds <- resource_discrete(3)
#' 
#' # Check initial state
#' beds$size()      # 3
#' beds$n_free()    # 3
#' beds$queue_size() # 0
#' 
#' # Block resources
#' i <- 101; curtime <- 0.0
#' beds$attempt_block()  # Uses i and curtime from environment
#' 
#' # Or explicitly
#' beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
#' 
#' # Check patient status
#' beds$is_patient_using(101)     # TRUE
#' beds$is_patient_in_queue(102)  # FALSE
#' 
#' @export
resource_discrete <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || n < 0 || n != as.integer(n)) {
    stop("n must be a single integer >= 0")
  }
  
  # Create the environment
  env <- new.env()
  
  # Create the C++ object using XPtr
  env$.ptr <- create_discrete_resource_cpp(as.integer(n))
  
  # Size method
  env$size <- function() {
    discrete_resource_size_cpp(env$.ptr)
  }
  
  # Queue size method
  env$queue_size <- function() {
    discrete_resource_queue_size_cpp(env$.ptr)
  }
  
  # Number of free resources method
  env$n_free <- function() {
    discrete_resource_n_free_cpp(env$.ptr)
  }
  
  # Get patients using resource
  env$patients_using <- function() {
    discrete_resource_patients_using_cpp(env$.ptr)
  }
  
  # Get start times of patients using resource
  env$patients_using_times <- function() {
    discrete_resource_patients_using_times_cpp(env$.ptr)
  }
  
  # Get queue start times
  env$queue_start_times <- function() {
    discrete_resource_queue_start_times_cpp(env$.ptr)
  }
  
  # Get queue priorities
  env$queue_priorities <- function() {
    discrete_resource_queue_priorities_cpp(env$.ptr)
  }
  
  # Get full queue information as data.frame
  env$queue_info <- function(n = NULL) {
    if (is.null(n)) n <- env$queue_size()
    if (n <= 0) return(data.frame(patient_id = integer(0), priority = integer(0), start_time = numeric(0)))
    
    patient_ids <- env$next_patient_in_line(n)
    priorities <- env$queue_priorities()[1:length(patient_ids)]
    start_times <- env$queue_start_times()[1:length(patient_ids)]
    
    data.frame(
      patient_id = patient_ids,
      priority = priorities,
      start_time = start_times,
      stringsAsFactors = FALSE
    )
  }
  
  # Check if patient is in queue
  env$is_patient_in_queue <- function(patient_id) {
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    discrete_resource_is_patient_in_queue_cpp(env$.ptr, as.integer(patient_id))
  }
  
  # Check if patient is using resource
  env$is_patient_using <- function(patient_id) {
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    discrete_resource_is_patient_using_cpp(env$.ptr, as.integer(patient_id))
  }
  
  # Attempt to block a resource unit
  env$attempt_block <- function(patient_id = NULL, priority = 1L, start_time = NULL) {
    # Get patient_id from parent frame if not provided
    if (is.null(patient_id)) {
      patient_id <-get("i", envir = parent.frame(), inherits = TRUE)
    }
    
    # Get start_time from parent frame if not provided
    if (is.null(start_time)) {
      start_time <- get("curtime", envir = parent.frame(), inherits = TRUE)
    }
    
    # Validate inputs
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    if (!is.numeric(priority) || length(priority) != 1) {
      stop("priority must be a single number")
    }
    if (!is.numeric(start_time) || length(start_time) != 1) {
      stop("start_time must be a single number")
    }
    
    discrete_resource_attempt_block_cpp(env$.ptr, as.integer(patient_id), as.integer(priority), as.numeric(start_time))
  }
  
  # Free a resource unit
  env$attempt_free <- function(patient_id = NULL, remove_all = FALSE) {
    # Get patient_id from parent frame if not provided
    if (is.null(patient_id)) {
      patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
    }
    
    # Validate inputs
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    if (!is.logical(remove_all) || length(remove_all) != 1) {
      stop("remove_all must be a single logical value")
    }
    
    discrete_resource_attempt_free_cpp(env$.ptr, as.integer(patient_id), remove_all)
    invisible(NULL)
  }
  
  # Free a resource unit only if patient is using it
  env$attempt_free_if_using <- function(patient_id = NULL, remove_all = FALSE) {
    # Get patient_id from parent frame if not provided
    if (is.null(patient_id)) {
      patient_id <- get("i", envir = parent.frame(), inherits = TRUE)
    }
    
    # Validate inputs
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    if (!is.logical(remove_all) || length(remove_all) != 1) {
      stop("remove_all must be a single logical value")
    }
    
    discrete_resource_attempt_free_if_using_cpp(env$.ptr, as.integer(patient_id), remove_all)
    invisible(NULL)
  }
  
  # Get next patients in line
  env$next_patient_in_line <- function(n = 1L) {
    if (!is.numeric(n) || length(n) != 1 || n < 1) {
      stop("n must be a single positive integer")
    }
    
    discrete_resource_next_patient_in_line_cpp(env$.ptr, as.integer(n))
  }
  
  # Modify priority of a patient in queue
  env$modify_priority <- function(patient_id, new_priority) {
    # Validate inputs
    if (!is.numeric(patient_id) || length(patient_id) != 1) {
      stop("patient_id must be a single number")
    }
    if (!is.numeric(new_priority) || length(new_priority) != 1) {
      stop("new_priority must be a single number")
    }
    
    discrete_resource_modify_priority_cpp(env$.ptr, as.integer(patient_id), as.integer(new_priority))
    invisible(NULL)
  }
  
  # Add resource units
  env$add_resource <- function(n_to_add) {
    if (!is.numeric(n_to_add) || length(n_to_add) != 1 || n_to_add < 1) {
      stop("n_to_add must be a single positive integer")
    }
    
    discrete_resource_add_resource_cpp(env$.ptr, as.integer(n_to_add))
    invisible(NULL)
  }
  
  # Remove resource units
  env$remove_resource <- function(n_to_remove, current_time = NULL) {
    # Get current_time from parent frame if not provided
    if (is.null(current_time)) {
      current_time <- get("curtime", envir = parent.frame(), inherits = TRUE)
    }
    
    if (!is.numeric(n_to_remove) || length(n_to_remove) != 1 || n_to_remove < 1) {
      stop("n_to_remove must be a single positive integer")
    }
    
    if (!is.numeric(current_time) || length(current_time) != 1) {
      stop("current_time must be a single number")
    }
    
    discrete_resource_remove_resource_cpp(env$.ptr, as.integer(n_to_remove), as.numeric(current_time))
    invisible(NULL)
  }
  
  # Set class
  class(env) <- "resource_discrete"
  
  return(env)
}

#' Print method for resource_discrete
#' @param x A resource_discrete object
#' @param ... Additional arguments (ignored)
#' 
#' @keywords internal
#' 
#' @export
print.resource_discrete <- function(x, ...) {
  cat("Discrete Resource:\n")
  cat("  Total capacity:", x$size(), "\n")
  cat("  Free units:", x$n_free(), "\n")
  cat("  Queue size:", x$queue_size(), "\n")
  cat("  Patients using:", length(x$patients_using()), "\n")
  invisible(x)
}


# Shared input ------------------------------------------------------------

#' Shared input object
#'
#' Constructor for a lightweight "shared or immutable" value holder.
#'
#' `shared_input()` produces a simple object that wraps a value with
#' controlled mutability semantics. It can operate in two distinct modes:
#'
#' - **Immutable (non-shared)**: every modification produces a fresh, independent
#'   copy of the object (safe for parallel or functional code).
#' - **Shared (constrained)**: the object’s value is stored in a common
#'   environment shared across all aliases (by-reference semantics). This
#'   allows coordinated updates across multiple handles.
#'
#' The mode is determined either by the explicit argument `constrained`, or
#' by inheriting the value of a `constrained` variable in the parent frame.
#'
#' @param expr A value or expression to initialize the shared input with.
#'   The expression is evaluated immediately.
#' @param constrained Logical. If `TRUE`, creates a shared environment-backed
#'   object. If `FALSE`, creates an immutable copy-on-modify object.
#'   If `NULL` (default), the function looks up `constrained` in the calling
#'   environment; only an explicit `TRUE` enables shared mode.
#'
#' @return An object of class `shared_input_val` (immutable mode) or
#'   `shared_input_env` (shared mode), both inheriting from class
#'   `"shared_input"`. Each instance exposes the following user methods:
#'
#' \describe{
#'   \item{$value()}{Returns the current stored value.}
#'   \item{$modify(new_v)}{
#'     In immutable mode: returns a new independent wrapper with updated value.  
#'     In shared mode: updates the shared value by reference and returns a new
#'     wrapper pointing to the same shared state.
#'   }
#'   \item{$clone()}{Returns a deep copy (independent wrapper and independent
#'     internal state). Subsequent modifications on clones do not affect the
#'     original object or its aliases.}
#'   \item{$reset()}{Returns a new wrapper whose value is restored to the
#'     original initialization value. In both modes this creates an independent
#'     fresh state.}
#'   \item{$fork(n)}{Creates `n` independent deep clones as a list. Useful for
#'     generating multiple isolated copies quickly.}
#' }
#'
#' @details
#' - In **immutable mode**, each wrapper stores its value in closures
#'   (`make_val()`) and is fully copy-on-modify. No references are shared.
#' - In **shared mode**, all wrappers produced by `$modify()` or direct aliasing
#'   point to the same underlying environment (`state`). This means updating one
#'   updates all aliases until a `$clone()` or `$reset()` breaks the link.
#'
#' The underlying `state` environments are internal. Users should rely only on
#' the public methods above.
#'
#' Note: if the stored value itself is a reference type (e.g., environment,
#' external pointer, R6 object), those internal references remain shared
#' regardless of mode, following normal R semantics.
#'
#' @examples
#' # --- Immutable (default) mode ---
#' a <- shared_input(5)
#' a$value()                 # 5
#' a2 <- a$modify(a$value() + 7)
#' a$value()                 # 5
#' a2$value()                # 12
#' 
#' # Cloning and resetting
#' a3 <- a2$clone()
#' a4 <- a2$reset()
#' a3$value(); a4$value()    # 12, 5
#' 
#' # Forking
#' forks <- a$fork(3)
#' vapply(forks, function(x) x$value(), numeric(1))
#'
#' # --- Shared (constrained) mode ---
#' constrained <- TRUE
#' b1 <- shared_input(10)
#' b2 <- b1        # alias (same state)
#' b1$modify(11)
#' b1$value(); b2$value()  # both 11
#'
#' b3 <- b1$clone()
#' b1$modify(99)
#' b1$value(); b3$value()  # 99, 11
#'
#' # Reset breaks sharing
#' b4 <- b1$reset()
#' b4$value()              # 10
#' 
#'@export
shared_input <- function(expr, constrained = NULL) {
  # Resolve constrained default from parent env, but only TRUE counts as TRUE
  if (is.null(constrained)) {
    constrained_parent <- get0("constrained",
                               envir = parent.frame(),
                               inherits = TRUE,
                               ifnotfound = NA)
    constrained <- isTRUE(constrained_parent)
  } else {
    constrained <- isTRUE(constrained)
  }
  
  val <- expr  # already evaluated before call
  
  # ---------- immutable (non-shared)
  if (!constrained) {
    make_val <- function(v, init = v) {
      force(v); force(init)
      self <- list()
      self$value  <- local({ vv <- v; function() vv })
      self$modify <- function(new_v) make_val(new_v, init = init)
      self$clone  <- function()     make_val(v, init = init)  # deep copy
      self$reset  <- function()     make_val(init, init = init)
      self$fork   <- function(n)    replicate(n, self$clone(), simplify = FALSE)
      class(self) <- c("shared_input_val", "shared_input")
      self
    }
    return(make_val(val, init = val))
  }
  
  # ---------- shared (environment-backed with shared state)
  # All wrappers to the same object share `state`; clone breaks sharing.
  make_env <- local({
    make <- function(state) {
      self <- new.env(parent = emptyenv())
      self$`.__state__` <- state
      self$value  <- function() state$val
      # Functional modify: updates shared state, returns a fresh wrapper
      self$modify <- function(new_v) { state$val <- new_v; make(state) }
      # Deep clone: new independent state
      self$clone  <- function() {
        new_state <- new.env(parent = emptyenv())
        new_state$val  <- state$val
        new_state$init <- state$init
        make(new_state)
      }
      # Reset to initial (independent fresh state)
      self$reset  <- function() {
        new_state <- new.env(parent = emptyenv())
        new_state$val  <- state$init
        new_state$init <- state$init
        make(new_state)
      }
      # Make n independent deep clones quickly
      self$fork   <- function(n) replicate(n, self$clone(), simplify = FALSE)
      
      class(self) <- c("shared_input_env", "shared_input")
      self
    }
    make
  })
  
  state <- new.env(parent = emptyenv())
  state$val  <- val
  state$init <- val
  make_env(state)
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

  evt_r <- list(list(react=input_sub))
  
  names(evt_r) <- paste(name_evt)

  if (is.null(data_list)) {
    data_list <- evt_r
  } else{
    data_list <- append(data_list,evt_r)
  }


  return(data_list)
}

# Random stream of uniform numbers -------------------------------------------------------------------------------------------------------------------------------------------

#' Creates an environment (similar to R6 class) of random uniform numbers to be drawn from
#'
#' @param stream_size Length of the vector of random uniform values to initialize
#'
#' @return Self (environment) behaving similar to R6 class
#'
#' @export
#'
#' @details
#' This function creates an environment object that behaves similar to an R6 class
#' but offers more speed vs. an R6 class.
#' 
#' The object is always initialized (see example below) to a specific vector of
#'  random uniform values. The user can then call the object with `obj$draw_number(n)`,
#'  where n is an integer, and will return the first n elements of the created
#'  vector of uniform values. It will automatically remove those indexes from the
#'  vector, so the next time the user calls `obj$draw_n()` it will already consider
#'  the next index.
#'  
#' The user can also access the latest elements drawn by accessing `obj$random_n` 
#'  (useful for when performing a luck adjustment), the current stream still 
#'  to be drawn using `obj$stream` and the original size (when created) using 
#'  `obj$stream_size`.
#'  
#'  If performing luck adjustment, the user can always modify the random value
#'  by using `obj$random_n <- luck_adj(...)` (only valid if used with the expression
#'  approach, not with `modify_item`)
#' 
#' 
#' @examples
#' stream_1 <- random_stream(1000)
#' number_1 <- stream_1$draw_n() #extract 1st index from the vector created
#' identical(number_1,stream_1$random_n) #same value
#' number_2 <- stream_1$draw_n() #gets 1st index (considers previous)
#' identical(number_2,stream_1$random_n) #same value

random_stream <- function(stream_size = 100) {
  self <- environment()
  
  # Initialize the stream with random numbers
  self$stream <- runif(stream_size)
  self$stream_size <- stream_size
  self$random_n <- numeric()
  
  # Function to generate a new stream of random numbers
  self$generate_stream <- function(size = self$stream_size) {
    self$stream <- runif(size)
  }
  # Function to draw a specified number of values from the stream and remove them
  self$draw_n <- function(n=1) {
    if (length(self$stream) < n) {
      warning("Stream is smaller than the number of numbers drawn. Generating a new stream of the correct size.")
      self$generate_stream(n)
      self$stream_size <- n
    }
    seq_index <- seq_len(n)
    self$random_n <- self$stream[seq_index]
    self$stream <- self$stream[-seq_index]
    return(self$random_n)
  }
  
  return(self)
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



#' Adjusted Value Calculation
#'
#' This function calculates an adjusted value over a time interval with optional discounting.
#' This is useful for instances when adding cycles may not be desirable, so one can perform
#' "cycle-like" calculations without needing cycles, offering performance speeds. See
#' the vignette on avoiding cycles for an example in a model.
#'
#' @param curtime Numeric. The current time point.
#' @param nexttime Numeric. The next time point. Must be greater than or equal to `curtime`.
#' @param by Numeric. The step size for evaluation within the interval.
#' @param expression An expression evaluated at each step. Use `.time` as the variable within the expression.
#' @param discount Numeric or NULL. The discount rate to apply, or NULL for no discounting.
#' @param vectorized_f boolean, FALSE by default. If TRUE, evaluates the expression once using `.time` as a vector.
#' If FALSE, it repeatedly evaluates the expression with time as a single value (slower).
#' 
#' @export
#'
#' @return Numeric. The calculated adjusted value.
#' 
#' @importFrom stats weighted.mean
#' @importFrom utils head
#' 
#' @details
#'  The user can use the `.time` variable to select the corresponding time of the sequence being evaluated.
#'  For example, in `curtime = 0, nexttime = 4, by = 1`, `.time` would correspond to `0, 1, 2, 3`.
#'  If using `nexttime = 4.2`, `0, 1, 2, 3, 4`
#' 
#'
#' @examples
#' # Define a function or vector to evaluate
#' bs_age <- 1
#' vec <- 1:8/10
#' 
#' # Calculate adjusted value without discounting
#' adj_val(0, 4, by = 1, expression = vec[floor(.time + bs_age)])
#' adj_val(0, 4, by = 1, expression = .time * 1.1)
#' #same result since .time * 1.1 can be vectorized w.r.t time
#' adj_val(0, 4, by = 1, expression = .time * 1.1, vectorized_f = TRUE)
#'
#' # Calculate adjusted value with discounting
#' adj_val(0, 4, by = 1, expression = vec[floor(.time + bs_age)], discount = 0.03)
#' 
adj_val <- function(curtime, nexttime, by, expression, discount = NULL, vectorized_f = FALSE) {
  duration <- nexttime - curtime
  if (duration < 0) stop("curtime - nexttime is negative (negative duration)")
  if (duration == 0) return(0)
  if (!is.null(discount) && is.infinite(discount)) return(0)
  
  n_steps <- floor(duration / by)
  times <- curtime + by * seq(0, n_steps)
  if (length(times) == 0 || tail(times, 1) < nexttime) {
    times <- c(times, nexttime)
  }
  
  intervals <- diff(times)
  eval_times <- head(times, -1)
  
  parent <- parent.frame()
  expr_sub <- substitute(expression)
  
  if(vectorized_f){
    values <- eval(expr_sub, envir = list(.time = eval_times), enclos = parent)
  } else{
    values <- vapply(eval_times, function(tt) {
      val <- eval(expr_sub, envir = list(.time = tt), enclos = parent)
      if (is.na(val)) stop("NA value encountered during evaluation at time: ", tt)
      val
    }, numeric(1))
  }
  
  if (!is.null(discount)) {
    # Compute present value of $1 over each interval (weights)
    weights <- disc_ongoing_v(
      lcldr = discount,
      lclprvtime = eval_times,
      lclcurtime = times[-1],
      lclval = rep(1, length(eval_times))
    )
    return(weighted.mean(values, w = weights))
  } else {
    return(sum(values * intervals) / duration)
  }
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
#' @importFrom lifecycle deprecate_warn
#'
#' @examples 
#' disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' 
#'
#' @export

disc_ongoing <- function(lcldr=0.035, lclprvtime, lclcurtime, lclval){
  lifecycle::deprecate_warn("2.0.0","disc_ongoing()","disc_ongoing_v()")
  
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



#' Calculate instantaneous discounted costs or qalys
#'
#' @param lcldr The discount rate
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#' 
#' @importFrom lifecycle deprecate_warn
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
  lifecycle::deprecate_warn("2.0.0","disc_instant()","disc_instant_v()")
  
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
#' @importFrom lifecycle deprecate_warn
#' 
#' @details
#' Note this function counts both extremes of the interval, so the example below
#'  would consider 25 cycles, while disc_cycle_v leave the right interval open
#'
#'
#' @examples 
#' disc_cycle(lcldr=0.035, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' 
#'
#' @export

disc_cycle <- function(lcldr=0.035, lclprvtime=0, cyclelength,lclcurtime, lclval,starttime=0){
  lifecycle::deprecate_warn("2.0.0","disc_cycle()","disc_cycle_v()")
  
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
#'                  d <- c <- k <- 67
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
  # Keep atoms/symbols as-is
  if (is.null(ee) || is.atomic(ee) || is.symbol(ee) || is.name(ee)) return(ee)
  
  # Recurse into expressions and pairlists
  if (is.expression(ee) || is.pairlist(ee)) {
    return(lapply(as.list(ee), ast_as_list))
  }
  
  # Calls → convert to lists and recurse over children (incl. function position)
  if (is.call(ee)) {
    return(lapply(as.list(ee), ast_as_list))
  }
  
  # Fallback
  ee
}


#' Extracts items and events by looking into assignments, modify_event and new_event
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
#'   a <- sum(5+7)
#'   
#'   ggplot()
#'   
#'   data.frame(x=1,b=2)
#'   
#'   list(b=5)
#'   
#'   a <- list(s=7)
#'   
#'   
#'   j <- 6
#'   if(TRUE){modify_event(list(j=5))}
#'   
#'   l <- 9
#'   
#'   afsa=ifelse(TRUE,"asda",NULL)
#'   
#'   
#'   o_exn = o_exn + 1
#'   
#'   a = NULL
#'   
#'   b = if(a){"CZ"}else{"AW"}
#'   
#'   rnd_prob_exn_sev = runif(1)
#'   
#'   exn_sev = rnd_prob_exn_sev <= p_sev
#'   
#'   o_exn_mod = o_exn_mod + if(exn_sev) { 0 } else { 1 }
#'   
#'   o_exn_sev = o_exn_sev + if(exn_sev) { 1 } else { 0 }
#'   
#'   o_rec_time_without_exn = (o_exn == 0) * 1
#'   
#'   o_rec_time_without_exn_sev = (o_exn_sev == 0) * 1
#'   
#'   o_c_exn = if(exn_sev) { c_sev } else { c_mod }
#'   
#'   o_other_c_exn_mod = if(exn_sev) { 0 } else { c_mod }
#'   
#'   o_other_c_exn_sev = if(exn_sev) { c_sev } else { 0 }
#'   
#'   o_qloss_exn = -if(exn_sev) { q_sev } else { q_mod }
#'   
#'   o_other_qloss_exn_mod = -if(exn_sev) { 0 } else { q_mod }
#'   
#'   o_other_qloss_exn_sev = -if(exn_sev) { q_sev } else { 0 }
#'   
#'   o_qloss_cg_exn = -if(exn_sev) { q_cg_sev } else { q_cg_mod }
#'   
#'   o_other_qloss_cg_exn_mod = -if(exn_sev) { 0 } else { q_cg_mod }
#'   
#'   o_other_qloss_cg_exn_sev = -if(exn_sev) { q_cg_sev } else { 0 }
#'   
#'   o_q = utility
#'   
#'   o_other_q_gold1 = if(gold == 1) { utility } else { 0 }
#'   
#'   o_other_q_gold2 = if(gold == 2) { utility } else { 0 }
#'   
#'   o_other_q_gold3 = if(gold == 3) { utility } else { 0 }
#'   
#'   o_other_q_gold4 = if(gold == 4) { utility } else { 0 }
#'   
#'   o_other_q_on_dup = if(on_dup) { utility } else { 0 }
#'   
#'   n_exn = n_exn + 1
#'   
#'   n_exn_mod = n_exn_mod + (1 - exn_sev)
#'   
#'   n_exn_sev = n_exn_sev + exn_sev
#'   
#'   u_adj_exn_lt = u_adj_exn_lt + if(exn_sev) { u_adj_sev_lt } else { u_adj_mod_lt }
#'   
#'   utility = u_gold - u_adj_exn_lt - u_mace_lt
#'   
#'   o_rec_utility = utility
#'   
#'   rnd_exn = runif(1)
#'   
#'   
#'   if(a==1){
#'     a=list(6+b)
#'     
#'     modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#'   } else{
#'     modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
#'     if(a>6){
#'       a=8
#'     }
#'     
#'   }
#'   
#'   
#'   if (sel_resp_incl == 1 & on_dup == 1) {
#'     
#'     modify_event(list(e_response = curtime, z = 6))
#'     
#'   }
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
  
  # Only list-like nodes have interesting structure here
  if (is.list(node) || is.call(node)) {
    # If we got an actual call, convert to list-form so downstream is uniform
    node_list <- if (is.call(node)) lapply(as.list(node), ast_as_list) else node
    
    func_name <- .node_head_name(node_list)
    
    # Case 1: modify_* / new_event 
    if (.is_call_named(node_list, c("modify_item_seq", "modify_item", "modify_event", "new_event"))) {
      type <- if (.is_call_named(node_list, c("modify_item_seq", "modify_item"))) "item" else "event"
      
      # Usually the 2nd argument is the list of definitions
      list_expr <- node_list[[2L]]
      if (is.list(list_expr)) {
        definition <- unlist(extract_defs(list_expr), recursive = TRUE, use.names = TRUE)
        definition <- definition[!(names(definition) %in% c("", NA))]
        if (length(definition)) {
          results <- rbind(
            results,
            data.frame(
              name = names(definition),
              type = type,
              conditional_flag = conditional_flag,
              definition = unname(definition),
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
    
    # Case 2: assignment via <- or = (including env/R6 LHS like env$x, x[[k]], x@slot)
    if (.is_call_named(node_list, c("<-", "="))) {
      # For a <- b, the list is typically: list("<-", LHS, RHS)
      lhs <- node_list[[2L]]
      rhs <- node_list[[3L]]
      
      target_name <- .lhs_to_name(lhs)
      if (!is.na(target_name) && nzchar(target_name)) {
        # Rebuild RHS for readable definition
        rhs_expr <- expr_from_list(rhs)
        results <- rbind(
          results,
          data.frame(
            name = target_name,
            type = "item",
            conditional_flag = conditional_flag,
            definition = clean_output(rhs_expr),
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    # ----- Case 3: assignment via assign("var", value, envir = env) 
    if (.is_call_named(node_list, "assign")) {
      # Positionals: assign(name, value, ...)
      if (length(node_list) >= 3L) {
        varname_node <- node_list[[2L]]
        value_node   <- node_list[[3L]]
        
        varname <- if (is.character(varname_node) && length(varname_node) == 1L) {
          varname_node
        } else {
          .lhs_to_name(varname_node)  # allow assign(sym, ...)
        }
        
        if (is.character(varname) && length(varname) == 1L && nzchar(varname)) {
          results <- rbind(
            results,
            data.frame(
              name = varname,
              type = "item",
              conditional_flag = conditional_flag,
              definition = clean_output(expr_from_list(value_node)),
              stringsAsFactors = FALSE
            )
          )
        }
      }
      # We ignore envir=… for the output schema, but this no longer warns or errors.
    }
    
    # ----- Case 4: if() — mark children as conditional 
    if (.is_call_named(node_list, "if")) {
      conditional_flag <- TRUE
    }
    
    # ----- Recurse into ALL children, including function position
    for (child in node_list) {
      results <- rbind(results, extract_elements_from_list(child, conditional_flag))
    }
  }
  
  # Final tidy
  results <- results[!(is.na(results$name) | results$name == "" | is.na(results$definition)), , drop = FALSE]
  results
}


# Return a single operator/name for the "head" of a node, or NULL.
# Works for calls, lists produced by ast_as_list(), symbols, etc.
.node_head_name <- function(node) {
  if (is.null(node)) return(NULL)
  
  # If we still have a language object (call/symbol), handle directly
  if (is.call(node))       return(.node_head_name(node[[1L]]))
  if (is.symbol(node))     return(as.character(node))
  if (is.name(node))       return(as.character(node))
  
  # If this is the list-form AST (from ast_as_list)
  if (is.list(node) && length(node) >= 1L) {
    return(.node_head_name(node[[1L]]))
  }
  
  # Atoms: no head
  NULL
}

# Predicate to test if node is a call with a given name (or one of names)
.is_call_named <- function(node, nm) {
  fn <- .node_head_name(node)
  is.character(fn) && length(fn) == 1L && fn %in% nm
}

# Rebuild a language object from the list-form AST.
expr_from_list <- function(lst) {
  if (is.null(lst)) return(NULL)
  
  # Already a language object
  if (is.symbol(lst) || is.call(lst) || is.name(lst)) return(lst)
  
  # Atomic vector → keep as-is
  if (is.atomic(lst)) return(lst)
  
  # Pairlist/expression/lists → rebuild recursively
  if (is.list(lst)) {
    if (length(lst) == 0L) return(NULL)
    head_expr <- expr_from_list(lst[[1L]])
    args      <- lapply(lst[-1L], expr_from_list)
    # IMPORTANT: use as.call so the head may itself be a call (e.g., beds$n_free)
    return(as.call(c(list(head_expr), args)))
  }
  
  # Fallback
  lst
}

# Compact deparse + tidy up quotes/backticks/spaces for stable string output
clean_output <- function(called) {
  txt <- paste0(deparse(called, width.cutoff = 500L), collapse = "")
  txt <- gsub("`", "", txt, fixed = TRUE)
  txt <- gsub("\"", "'", txt, fixed = TRUE)
  txt <- gsub("    ", "", txt, fixed = TRUE)
  txt
}

# Convert a (possibly list-form) LHS node into a printable "target name".
# For plain symbols: "x"; for complex: "env$x", "x[[k]]", "x@slot", etc.
.lhs_to_name <- function(lhs) {
  if (is.null(lhs)) return(NA_character_)
  if (is.symbol(lhs) || is.call(lhs)) return(clean_output(lhs))
  if (is.list(lhs))  return(clean_output(expr_from_list(lhs)))
  if (is.character(lhs) && length(lhs) == 1L) return(lhs)
  NA_character_
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
  
  # Try to extract the base object name on the LHS (e.g., x from x$y, x[[...]], x[i], x@slot, slot(x, "slot"))
  base_symbol_from_lhs <- function(lhs) {
    if (is.symbol(lhs)) {
      return(as.character(lhs))
    }
    if (!is.call(lhs)) return(NULL)
    
    op <- lhs[[1L]]
    # Operators as symbols?
    if (is.symbol(op)) {
      op_chr <- as.character(op)
      
      # x$y, x[[...]], x[...], x@slot
      if (op_chr %in% c("$", "[[", "[", "@")) {
        base <- lhs[[2L]]
        if (is.symbol(base)) return(as.character(base))
        return(NULL)
      }
      
      # slot(x, "foo") <- ...
      if (op_chr == "slot" && length(lhs) >= 2L) {
        base <- lhs[[2L]]
        if (is.symbol(base)) return(as.character(base))
        return(NULL)
      }
      
      # pkg::name  or pkg:::name  (not assignable in R) → ignore as target but don't error
      if (op_chr %in% c("::", ":::")) {
        return(NULL)
      }
    }
    
    # Function-call LHS of other kinds → we don't treat as assignable base names
    NULL
  }
  
  walk <- function(node) {
    if (is.null(node)) return(invisible())
    
    # Walk pairlists or expressions
    if (is.pairlist(node) || is.expression(node)) {
      lapply(as.list(node), walk)
      return(invisible())
    }
    
    if (is.call(node)) {
      fn <- node[[1L]]
      
      # <- and =
      if (is.symbol(fn) && (identical(fn, as.name("<-")) || identical(fn, as.name("=")))) {
        lhs <- node[[2L]]
        
        # Record base symbol if we can determine it
        nm <- base_symbol_from_lhs(lhs)
        if (!is.null(nm)) assigned <<- c(assigned, nm)
        
        # Recurse both sides to catch nested assignments
        walk(lhs)
        if (length(node) >= 3L) walk(node[[3L]])
      }
      
      # assign("var", value)
      else if (is.symbol(fn) && identical(fn, as.name("assign"))) {
        if (length(node) >= 2L) {
          v <- node[[2L]]
          if (is.character(v) && length(v) == 1L) {
            assigned <<- c(assigned, v)
          } else if (is.symbol(v)) {
            assigned <<- c(assigned, as.character(v))
          }
        }
        if (length(node) >= 3L) walk(node[[3L]])
      }
      
      # If the "function position" itself is a call, walk it (e.g., beds$n_free(), pkg::fun(), pkg:::fun())
      if (is.call(fn)) walk(fn)
      
      # Walk all arguments
      if (length(node) > 1L) {
        for (i in 2L:length(node)) walk(node[[i]])
      }
      
      return(invisible())
    }
    
    invisible()
  }
  
  walk(expr)
  unique(assigned)
}