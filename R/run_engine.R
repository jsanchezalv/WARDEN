#' Run the simulation engine without using a parallel engine, facilitating debugging
#'
#' @param arm_list A vector of the names of the interventions evaluated in the simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param input_list A list of all other inputs: drc, drq, psa_bool, init_event_list, evt_react_list,
#' uc_lists = list(util_ongoing_list,util_instant_list,util_cycle_list,cost_ongoing_list,cost_instant_list,cost_cycle_list),
#' input_out,ipd,arm_list,simulation,npats,n_sim
#' @param pb progress bar
#' @param seed Starting seed to be used for the whole analysis
#'
#' @return A data frame with the simulation results
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom data.table rbindlist
#'
#' @noRd



run_engine <- function(arm_list,
                            common_pt_inputs=NULL,
                            unique_pt_inputs=NULL,
                            common_arm_inputs = NULL,
                            input_list = NULL,
                            pb = pb,
                            seed = seed){
  # Initial set-up --------------------------
  arm_list <- arm_list
  simulation <- input_list$simulation
  sens <- input_list$sens
  n_sensitivity <- input_list$n_sensitivity
  length_sensitivities <- input_list$length_sensitivities
  n_sim <- input_list$n_sim
  n_arms <- length(arm_list)
  npats <- input_list$npats
  psa_bool <- input_list$psa_bool
  env_setup_pt <- input_list$env_setup_pt
  env_setup_arm <- input_list$env_setup_arm
  env_setup_arm_common <- input_list$env_setup_arm_common
  

  #1 Loop per patient ----------------------------------------------------------
  patdata <- vector("list", length=npats) # empty list with npats elements

  temp_log_pt <- list()

    
  tryCatch({
    
    list_discrete_resources <- list()
    for (obj in ls(input_list)) {
      if(inherits(input_list[[obj]],"resource_discrete")){
        new_obj <- list(input_list[[obj]])
        names(new_obj) <- obj
        list_discrete_resources <- c(list_discrete_resources,new_obj)
      }
    }
    
    l_disres <- length(list_discrete_resources)
    if(l_disres>0){
      names_disres <- names(list_discrete_resources)
      cloned_resources <- list()
      for (obj in 1:l_disres) {
        cloned_resources[[obj]] <- discrete_resource_clone(list_discrete_resources[[obj]], n_arms)
      }
    }
  
  for (i in 1:npats) {
    set.seed((simulation*1007 + i*53) * seed)
    if((((sens - 1) * n_sim * npats) + ((simulation - 1) * npats) + i) %% ceiling(npats*n_sim*length_sensitivities / min(npats*length_sensitivities*n_sim,50)) == 0){
      pb(sprintf("Simulation %g", simulation))
      }
    
    #Create empty pat data for each arm
    this_patient <- list()
    input_list_pt <- rlang::env_clone(input_list , parent.env(input_list))
    input_list_pt$i <- i
    
    #Extract the inputs that are common for each patient across interventions
    if(!is.null(common_pt_inputs)){
      
      if(env_setup_pt){
        load_inputs2(inputs = input_list_pt,list_uneval_inputs = common_pt_inputs)
      } else{
        input_list_pt <- as.environment(
          load_inputs(inputs = as.list(input_list_pt),
                      list_uneval_inputs = common_pt_inputs)
          )
        parent.env(input_list_pt) <- parent.env(input_list)
        
      }
      
      
      if(input_list$debug){ 
        dump_info <- debug_inputs(input_list,input_list_pt)
        
        
        names(dump_info) <- paste0("Analysis: ", input_list_pt$sens," ", input_list_pt$sens_name_used,
                                   "; Sim: ", input_list_pt$simulation,
                                   "; Patient: ", input_list_pt$i,
                                   "; Initial Patient Conditions"
        )
        
        temp_log_pt <- c(temp_log_pt,dump_info)
      }
    }
    
    #2 Loop per treatment ------------------------------------------------------
    temp_log <- list()
    
    for (arm in arm_list) {
      
      #Extract the inputs that are unique for each patient-intervention
      input_list_arm <- rlang::env_clone(input_list_pt , parent.env(input_list_pt))
      input_list_arm$arm <- arm
      
      if(l_disres>0){
        which_arm <- which(arm==arm_list)
        for (obj in 1:l_disres) {
          input_list_arm[[names_disres[[obj]]]] <-  cloned_resources[[obj]][[which_arm]]
        }
      }
      
      #Unique patient-arm inputs
      set.seed(seed*(simulation*1007 + i*53 + which(arm==arm_list)))
      # set.seed(seed*(simulation*1007 + i*191))
      # Initialize values to prevent errors

      if(!is.null(unique_pt_inputs)){
        
        if(env_setup_arm){
          load_inputs2(inputs = input_list_arm,list_uneval_inputs = unique_pt_inputs)
        } else{
          input_list_arm <- as.environment(
            load_inputs(inputs = as.list(input_list_arm),
                        list_uneval_inputs = unique_pt_inputs)
            )
          parent.env(input_list_arm) <- parent.env(input_list_pt)
        }
        
        if(input_list_pt$debug){ 
          dump_info <- debug_inputs(input_list_pt,input_list_arm)
          
          names(dump_info) <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                                     "; Sim: ", input_list_arm$simulation,
                                     "; Patient: ", input_list_arm$i,
                                     "; Initial Patient-Arm Conditions"
          )
          
          temp_log <- c(temp_log,dump_info)
        }
      }
      
      # Generate event list
      set.seed(seed*(simulation*1007 + i*349))
      
      if (!is.null(input_list_arm$init_event_list)) {
        priority_order <- input_list_arm$init_event_list[[1]]$evts
      } else {
        priority_order <- "start"  # default if no events defined
      }
      event_queue <- queue_create(priority_order)
      
      if (is.null(input_list_arm$init_event_list)) {
        # No events defined, add start event at time 0
        start_events <- setNames(0, "start")
        new_event2(start_events, event_queue, i)
      } else {
        # Generate initial events
        evt_list <- do.call("initiate_evt", list(arm, input_list_arm))
        
        if (input_list_arm$debug) {
          names_input <- names(evt_list$time_data)
          prev_value <- setNames(vector("list", length(names_input)), names_input)
          prev_value[names_input] <- evt_list$time_data[names_input]
          prev_value["cur_evtlist"] <- list(setNames(rep(Inf, length(input_list_arm$init_event_list[[1]]$evts)), 
                                                     input_list_arm$init_event_list[[1]]$evts))
          dump_info <- list(
            list(
              prev_value = prev_value,
              cur_value = c(evt_list[["time_data"]], evt_list["cur_evtlist"])
            )
          )
          names(dump_info) <- paste0("Analysis: ", input_list_arm$sens, " ", input_list_arm$sens_name_used,
                                     "; Sim: ", input_list_arm$simulation,
                                     "; Patient: ", input_list_arm$i,
                                     "; Initialize Time to Events for Patient-Arm")
          temp_log <- c(temp_log, dump_info)
        }
        
        # Load time data into patient-arm environment
        list2env(as.list(evt_list$time_data), input_list_arm)
        
        # Add events to the shared event queue
        if (length(evt_list$cur_evtlist) > 0) {
          new_event2(evt_list$cur_evtlist, event_queue, i)
        }
      }

      # 3 Loop per event --------------------------------------------------------

      this_patient[[arm]]$evtlist <- NULL

      input_list_arm$curtime <- 0
      
      if(input_list$accum_backwards){
        input_list_arm$ongoing_inputs_lu <- paste0(input_list_arm$uc_lists$ongoing_inputs,"_lastupdate",recycle0 = TRUE)
        input_out_v <- c(input_list_arm$input_out,
                         input_list_arm$ongoing_inputs_lu
        )
      }else{
        input_out_v <- c(input_list_arm$input_out)
      }

      n_evt <- 0
   
      if(input_list$accum_backwards){
        input_list_arm$ongoing_inputs_lu <- paste0(input_list_arm$uc_lists$ongoing_inputs,"_lastupdate",recycle0 = TRUE)
        inputs_out_v <- c(input_list_arm$input_out,
                          input_list_arm$ongoing_inputs_lu
        )
      }else{
        inputs_out_v <-  input_list_arm$input_out
      }

      # Update the event queue reference for new_event2, modify_event2, etc.
      # Note that we are only assigning the pointer, so any changes to input_list_arm$cur_evtlist
      # will equally affect event_queue
      assign("cur_evtlist", event_queue, envir = input_list_arm)

      while (!queue_empty(event_queue)) {
        if(is.infinite(next_event(1,event_queue)$time)){
          break
        }
        # Get next event
        next_evt <- pop_and_return_event(event_queue)
        current_event <- next_evt$event_name
        current_time <- next_evt$time
        
        # Calculate prevtime correctly (current curtime becomes prevtime)
        current_prevtime <- input_list_arm$curtime
        if(is.infinite(current_prevtime)){next}
        
        n_evt <- n_evt + 1
        
        input_list_arm <- react_evt(list(evt = current_event,
                                         curtime = current_prevtime,
                                         evttime = current_time),
                                    arm,
                                    input_list_arm)
        
        #Get extra objects to be exported
        if(is.null(inputs_out_v)){
          extra_data <- list()
        } else{
          extra_data <-  mget(inputs_out_v, input_list_arm) 
        }
        extra_data <- extra_data[!vapply(extra_data, is.null, TRUE)]
          
              this_patient[[arm]]$evtlist[[n_evt]] <- c(evtname = current_event,
                                                        evttime = current_time,
                                                        pat_id = i,
                                                        arm = arm,
                                                        extra_data
              )
              

        }
      
      
      temp_log <- c(temp_log,input_list_arm$log_list)
    }
    temp_log_pt <- c(temp_log_pt,temp_log)

    patdata[[i]] <- this_patient
  }
  
  input_list$log_list <- lapply(temp_log_pt,transform_debug)
  
  
# Compute outputs ---------------------------------------------------------


  #Compute the outputs and format the data
  final_output <- compute_outputs(patdata, input_list)
  
  if(input_list$debug){
    final_output$log_list <- input_list$log_list
  }
    return(final_output)
  
  }, error = function(e) {

    if(input_list$debug){
      
      if(!is.null(input_list_arm$log_list)){
        temp_log <- c(temp_log,input_list_arm$log_list)
      }
      
      if(!is.null(temp_log)){
        temp_log_pt <- c(temp_log_pt,temp_log)
      }
      
      final_output <- list()  
      final_output$log_list <- lapply(temp_log_pt,transform_debug)
      final_output$error_m <- paste0(e$message," in ", e$call,
                                     paste0(". Error in patient: ", if(exists("i")){i},
                                             "; arm: ", if(exists("arm")){arm},
                                             "; event: ", if(exists("Evt")){Evt$evt},
                                             "; time: ", if(exists("Evt")){Evt$evttime})
                                     )
      return(final_output)
    }else if(input_list$continue_on_error){
      final_output <- list()  
      
      final_output$error_m <- paste0(e$message," in ", e$call,
                                     paste0(". Error in patient: ", if(exists("i")){i},
                                             "; arm: ", if(exists("arm")){arm},
                                             "; event: ", if(exists("Evt")){Evt$evt},
                                             "; time: ", if(exists("Evt")){Evt$evttime})
      )
      return(final_output)
      
    } else{
      
      stop(paste0(" ", e$message," in ", e$call,
                  paste0(". Error in patient: ", if(exists("i")){i},
                           "; arm: ", if(exists("arm")){arm},
                           "; event: ", if(exists("Evt")){Evt$evt},
                           "; time: ", if(exists("Evt")){Evt$evttime})
      ))
    }
  } )

}
