#' Run the simulation engine with constrained event processing using event queues
#'
#' This function implements a discrete event simulation engine that processes events
#' in chronological order across all patients within each arm, using shared event queues.
#' Unlike run_engine which processes patients sequentially, this function loads all
#' patients for an arm and then processes events in time order regardless of patient.
#'
#' @param arm_list A vector of the names of the interventions evaluated in the simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention.
#'   These inputs are loaded into patient-specific environments that are children of the arm-level environment.
#' @param unique_pt_inputs A list of inputs that change across each intervention.
#'   These inputs are loaded into patient-arm specific environments.
#' @param input_list A list of all other inputs including:
#'   \itemize{
#'     \item drc, drq: discount rates for costs and QALYs
#'     \item psa_bool: whether PSA is enabled
#'     \item init_event_list: list of initial events and event times
#'     \item evt_react_list: list of event reactions
#'     \item uc_lists: list containing util_ongoing_list, util_instant_list, etc.
#'     \item input_out: vector of output variable names
#'     \item npats: number of patients
#'     \item simulation: simulation number
#'   }
#' @param pb Progress bar function
#' @param seed Starting seed to be used for the whole analysis
#'
#' @return A list containing simulation results with the same structure as run_engine
#' 
#' @details 
#' The function creates a hierarchical environment structure per arm:
#' \enumerate{
#'   \item Simulation-level environment (cloned per arm)
#'   \item Patient-level environments (children of arm simulation environment)
#'   \item Patient-arm environments (children of patient environment)
#' }
#' 
#' This allows patients within an arm to modify shared inputs (like hospital resources)
#' while maintaining patient-specific inputs (like age) in isolation.
#' 
#' Events are processed using a shared event queue per arm, ensuring chronological
#' ordering across all patients within the arm.
#'
#' @importFrom purrr map
#' @importFrom purrr map_dbl  
#' @importFrom data.table rbindlist
#' @importFrom rlang env_clone
#'
#' @noRd
run_engine_constrained <- function(arm_list,
                                   common_pt_inputs = NULL,
                                   unique_pt_inputs = NULL,
                                   input_list = NULL,
                                   pb = pb,
                                   seed = seed) {
  
  
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
  debug <- input_list$debug
  
  # Get priority order for event queue (from init_event_list)
  if (!is.null(input_list$init_event_list)) {
    priority_order <- input_list$init_event_list[[1]]$evts
  } else {
    priority_order <- "start"  # default if no events defined
  }
  
  
  # Storage for patient data per arm
  patdata <- vector("list", length=npats) # empty list with npats elements
  
  temp_log_pt <- list()
  temp_log <- list()
  
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
    
    # 1 Loop per arm ----------------------------------------------------------
    for (arm in arm_list) {
      # Clone simulation environment for this arm
      input_list_arm_base <- rlang::env_clone(input_list, parent.env(input_list)) 
      input_list_arm_base$arm <- arm
      
      if(l_disres>0){
        which_arm <- which(arm==arm_list)
        for (obj in 1:l_disres) {
          input_list_arm_base[[names_disres[[obj]]]] <-  cloned_resources[[obj]][[which_arm]]
        }
      }
      
      # Create event queue for this arm
      event_queue <- queue_create(priority_order)
      
      # Storage for patient environments and data
      patient_arm_environments <- vector("list", length = npats)
      
      
      # 2 Load all patients for this arm and initialize events ---------------
      for (i in 1:npats) {
        set.seed((simulation * 1007 + i * 53) * seed)
        
        # Progress update
        if ((((sens - 1) * n_sim * npats * n_arms) + ((simulation - 1) * npats * n_arms) + i) %% 
            ceiling(npats * n_arms * n_sim * length_sensitivities / min(npats * n_arms * length_sensitivities * n_sim, 50)) == 0) {
          pb(sprintf("Simulation %g", simulation))
        }
        
        # Create patient environment (child of arm base environment)
        input_list_pt <- rlang::env_clone(input_list_arm_base, parent.env(input_list_arm_base))
        input_list_pt$i <- i
        
        # Load common patient inputs if they exist
        if (!is.null(common_pt_inputs)) {
          if (env_setup_pt) {
            load_inputs2(inputs = input_list_pt, list_uneval_inputs = common_pt_inputs)
          } else {
            input_list_pt <- as.environment(
              load_inputs(inputs = as.list(input_list_pt),
                          list_uneval_inputs = common_pt_inputs)
            )
            parent.env(input_list_pt) <- parent.env(input_list_arm_base)
          }
          
          if (input_list$debug) {
            dump_info <- debug_inputs(input_list_arm_base, input_list_pt)
            names(dump_info) <- paste0("Analysis: ", input_list_pt$sens, " ", input_list_pt$sens_name_used,
                                       "; Sim: ", input_list_pt$simulation,
                                       "; Arm: ", input_list_pt$arm,
                                       "; Patient: ", input_list_pt$i,
                                       "; Initial Patient Conditions")
            temp_log_pt <- c(temp_log_pt, dump_info)
          }
        }
        
        
        set.seed(seed*(simulation*1007 + i*53 + which(arm==arm_list)))
        
        # Create patient-arm environment (child of patient environment)
        input_list_arm <- rlang::env_clone(input_list_pt, parent.env(input_list_pt))
        
        # Load unique patient-arm inputs
        if (!is.null(unique_pt_inputs)) {
          if (env_setup_arm) {
            load_inputs2(inputs = input_list_arm, list_uneval_inputs = unique_pt_inputs)
          } else {
            input_list_arm <- as.environment(
              load_inputs(inputs = as.list(input_list_arm),
                          list_uneval_inputs = unique_pt_inputs)
            )
            parent.env(input_list_arm) <- parent.env(input_list_pt)
          }
          
          if (input_list_pt$debug) {
            dump_info <- debug_inputs(input_list_pt, input_list_arm)
            names(dump_info) <- paste0("Analysis: ", input_list_arm$sens, " ", input_list_arm$sens_name_used,
                                       "; Sim: ", input_list_arm$simulation,
                                       "; Arm: ", input_list_arm$arm,
                                       "; Patient: ", input_list_arm$i,
                                       "; Initial Patient-Arm Conditions")
            temp_log_pt <- c(temp_log_pt, dump_info)
          }
        }
        
        
        # Initialize events for this patient-arm
        set.seed(seed * (simulation * 1007 + i * 349))
        
        if (is.null(input_list_arm$init_event_list)) {
          # No events defined, add start event at time 0
          start_events <- setNames(0, "start")
          new_event(start_events, event_queue, i)
        } else {
          # Generate initial events
          evt_list <- do.call("initiate_evt", list(arm, input_list_arm))
          
          if (input_list_pt$debug) {
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
                                       "; Arm: ", arm,
                                       "; Initialize Time to Events for Patient-Arm")
            temp_log_pt <- c(temp_log_pt, dump_info)
          }
          
          # Load time data into patient-arm environment
          list2env(as.list(evt_list$time_data), input_list_arm)
          
          # Add events to the shared event queue
          if (length(evt_list$cur_evtlist) > 0) {
            new_event(evt_list$cur_evtlist, event_queue, i)
          }
        }
        
        # Initialize output list for this patient
        input_list_arm$curtime <- 0
        
        # Set up accumulators if using backwards accumulation
        if (input_list$accum_backwards) {
          input_list_arm$ongoing_inputs_lu <- paste0(input_list_arm$uc_lists$ongoing_inputs, "_lastupdate", recycle0 = TRUE)
          inputs_out_v <- c(input_list_arm$input_out, input_list_arm$ongoing_inputs_lu)
          
          # Initialize ongoing_list_temp for backwards accumulation
          if (!is.null(input_list_arm$uc_lists$ongoing_inputs)) {
            input_list_arm$ongoing_list_temp <- setNames(
              rep(0, length(input_list_arm$uc_lists$ongoing_inputs)),
              input_list_arm$uc_lists$ongoing_inputs
            )
          }
        } else {
          inputs_out_v <- c(input_list_arm$input_out)
        }
        
        input_list_arm$n_evt <- 0
        # Store patient-arm environment
        patient_arm_environments[[i]] <- input_list_arm
      }
      
      # 3 Process all events for this arm in chronological order -------------
      
      n_evt <- 0
      
      # Process events while queue is not empty
      while (!queue_empty(event_queue)) {
        
        if(is.infinite(next_event(1,event_queue)$time)){
          break
        }
        
        # Get next event
        next_evt <- pop_and_return_event(event_queue)
        current_patient_id <- next_evt$patient_id
        current_event <- next_evt$event_name
        current_time <- next_evt$time
        # Get the appropriate patient-arm environment
        input_list_arm <- patient_arm_environments[[current_patient_id]]

        # Calculate prevtime correctly (current curtime becomes prevtime)
        current_prevtime <- input_list_arm$curtime
        if(is.infinite(current_prevtime)){
          if(!queue_empty(event_queue)){next}else{break}
        }
        
        # Update the event queue reference for new_event, modify_event, etc.
        # Note that we are only assigning the pointer, so any changes to input_list_arm$cur_evtlist
        # will equally affect event_queue
        assign("cur_evtlist", event_queue, envir = input_list_arm)
        
        input_list_arm$n_evt <- input_list_arm$n_evt + 1
        
        input_list_arm <- react_evt(list(evt = current_event,
                                         curtime = current_prevtime,
                                         evttime = current_time),
                                    arm,
                                    input_list_arm)
       
        patient_arm_environments[[current_patient_id]] <-  input_list_arm 
        
        
        #Get extra objects to be exported
        if(is.null(inputs_out_v)){
          extra_data <- list()
        } else{
          extra_data <-  mget(inputs_out_v, input_list_arm) 
        }
        extra_data <- extra_data[!vapply(extra_data, is.null, TRUE)]
        
        
        patdata[[current_patient_id]][[arm]]$evtlist[[input_list_arm$n_evt]] <- c(evtname = current_event,
                                                  evttime = current_time,
                                                  pat_id = current_patient_id,
                                                  arm = arm,
                                                  extra_data
        )
       
        temp_log_pt <- c(temp_log_pt,input_list_arm$log_list)
        
      }
    }
    
    # Return patdata in the same structure as run_engine
    final_output <- compute_outputs(patdata, input_list)
    
    input_list$log_list <- lapply(temp_log_pt,transform_debug)
    
    # Add debugging information if enabled
    if(input_list$debug){
      final_output$log_list <- input_list$log_list 
    }
   
    
    return(final_output)
    
  }, error = function(e) {
    if(input_list$debug){
      
      final_output <- list()  
      final_output$log_list <- lapply(temp_log_pt,transform_debug)
      final_output$error_m <- paste0(e$message, " in ", paste0(e$call,collapse=", "),
                                     paste0(". Error in patient: ", if (exists("current_patient_id")) current_patient_id,
                                            "; arm: ", if (exists("arm")) arm,
                                            "; event: ", if (exists("current_event")) current_event,
                                            "; time: ", if (exists("current_time")) current_time))
      return(final_output)
    }else if(input_list$continue_on_error){
      final_output <- list()  
      
      final_output$error_m <- paste0(e$message, " in ", paste0(e$call,collapse=", "),
                                     paste0(". Error in patient: ", if (exists("current_patient_id")) current_patient_id,
                                            "; arm: ", if (exists("arm")) arm,
                                            "; event: ", if (exists("current_event")) current_event,
                                            "; time: ", if (exists("current_time")) current_time))
      return(final_output)
      
    } else{
      
      stop(paste0(e$message, " in ", paste0(e$call,collapse=", "),
                  paste0(". Error in patient: ", if (exists("current_patient_id")) current_patient_id,
                         "; arm: ", if (exists("arm")) arm,
                         "; event: ", if (exists("current_event")) current_event,
                         "; time: ", if (exists("current_time")) current_time))
      )
    }
  } )
}