#' Run the simulation engine without using a parallel engine, facilitating debugging
#'
#' @param arm_list A vector of the names of the interventions evaluated in the simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param input_list A list of all other inputs: drc, drq, psa_bool, init_event_list, evt_react_list,
#' uc_lists = list(util_ongoing_list,util_instant_list,util_cycle_list,cost_ongoing_list,cost_instant_list,cost_cycle_list),
#' input_out,ipd,arm_list,simulation,npats,n_sim
#' @param pb progress bar
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
                            input_list = NULL,
                       pb = pb){
  # Initial set-up --------------------------
  arm_list <- arm_list
  simulation <- input_list$simulation
  sens <- input_list$sens
  n_sensitivity <- input_list$n_sensitivity
  n_sim <- input_list$n_sim
  npats <- input_list$npats
  psa_bool <- input_list$psa_bool


  #1 Loop per patient ----------------------------------------------------------
  patdata <- vector("list", length=npats) # empty list with npats elements

  temp_log_pt <- list()

    
  tryCatch({
  
  for (i in 1:npats) {
    set.seed(sens*100037 + simulation*1007 + i*53)
    if((((sens - 1) * n_sim * npats) + ((simulation - 1) * npats) + i) %% ceiling(npats*n_sim*n_sensitivity / 50) == 0){
      pb(sprintf("Simulation %g", simulation))
      }
    
    #Create empty pat data for each arm
    this_patient <- list()
    input_list_pt <- c(input_list,list(i=i))
    
    #Extract the inputs that are common for each patient across interventions
    if(!is.null(common_pt_inputs)){
      for (inp in 1:length(common_pt_inputs)) {
        list.common_pt_inputs <- lapply(common_pt_inputs[inp],function(x) eval(x, input_list_pt))
        #If using pick_eval_v or other expressions, the lists are not deployed, so this is necessary to do so
        if(any(is.null(names(list.common_pt_inputs)), names(list.common_pt_inputs)=="") & length(list.common_pt_inputs)==1) {
          input_list_pt <- c(input_list_pt, list.common_pt_inputs[[1]])
        } else{
        if (!is.null(names(list.common_pt_inputs[[1]]))) {
          warning("Item ", names(list.common_pt_inputs), " is named. It is advised to assign unnamed objects if they are going to be processed in the model, as they can create errors depending on how they are used within the model.\n")
        }
        input_list_pt <- c(input_list_pt,list.common_pt_inputs)
        }
      }
      
      if(input_list$debug){ 
        names_pt_input <- names(common_pt_inputs)
        prev_value <- setNames(vector("list", length(common_pt_inputs)), names_pt_input)
        prev_value[names_pt_input] <- input_list[names_pt_input]
        dump_info <- list(
          list(
            prev_value = prev_value,
            cur_value  = input_list_pt[names_pt_input]
          )
        )
        
        names(dump_info) <- paste0("Analysis: ", input_list_pt$sens," ", input_list_pt$sens_name_used,
                                   "; Sim: ", input_list_pt$sim,
                                   "; Patient: ", input_list_pt$i,
                                   "; Initial Patient Conditions"
        )
        
        temp_log_pt <- c(temp_log_pt,dump_info)
      }
      
    }

    #Make sure there are no duplicated inputs in the model, if so, take the last one
    duplic <- duplicated(names(input_list_pt),fromLast = T)
    if (sum(duplic)>0 & i==1 & simulation==1 & sens==1) { warning("Duplicated items detected in the Patient, using the last one added.\n")  }
    input_list_pt <- input_list_pt[!duplic]

    #2 Loop per treatment ------------------------------------------------------
    temp_log <- list()
    
    for (arm in arm_list) {
      
      set.seed(sens*1037 + simulation*1007 + i*53 + which(arm==arm_list))
      # Initialize values to prevent errors
      output_list <- list(curtime = 0)
      
      #Extract the inputs that are unique for each patient-intervention
      input_list_arm <- NULL
      input_list_arm <- c(input_list_pt,list(arm=arm))
      
      
      if(!is.null(unique_pt_inputs)){
        for (inp in 1:length(unique_pt_inputs)) {
          list.unique_pt_inputs <- lapply(unique_pt_inputs[inp],function(x) eval(x, input_list_arm))
          #If using pick_eval_v or other expressions, the lists are not deployed, so this is necessary to do so
          if(any(is.null(names(list.unique_pt_inputs)), names(list.unique_pt_inputs)=="") & length(list.unique_pt_inputs)==1) {
            input_list_arm <- c(input_list_arm, list.unique_pt_inputs[[1]])
          } else{
          if ((!is.null(names(list.unique_pt_inputs[[1]]))) & i==1 & simulation==1 & sens==1) {
            warning("Item ", names(list.unique_pt_inputs), " is named. It is advised to assign unnamed objects if they are going to be processed in the model, as they can create errors depending on how they are used within the model.\n")
          }
          input_list_arm <- c(input_list_arm,list.unique_pt_inputs)
          }
        }
        
        if(input_list_pt$debug){ 
          names_pt_input <- names(unique_pt_inputs)
          prev_value <- setNames(vector("list", length(unique_pt_inputs)), names_pt_input)
          prev_value[names_pt_input] <- input_list_pt[names_pt_input]
          dump_info <- list(
            list(
              prev_value = prev_value,
              cur_value  = input_list_arm[names_pt_input]
            )
          )

          names(dump_info) <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                                     "; Sim: ", input_list_arm$sim,
                                     "; Patient: ", input_list_arm$i,
                                     "; Initial Patient-Arm Conditions"
          )
          
          temp_log <- c(temp_log,dump_info)
        }
      }

      #Make sure there are no duplicated inputs in the model, if so, take the last one
      duplic <- duplicated(names(input_list_arm),fromLast = T)
      if (sum(duplic)>0 & i==1 & simulation==1 & sens==1) { warning("Duplicated items detected in the Arm, using the last one added.\n")  }
      input_list_arm <- input_list_arm[!duplic]

      # Generate event list
      #if noeventlist, then just make start at 0
      if (is.null(input_list_arm$init_event_list)) {
        evt_list <- list(cur_evtlist = setNames(0,"start"), time_data = NULL)
      } else{
        evt_list <- do.call("initiate_evt",list(arm,input_list_arm))
      }

    if(input_list_pt$debug){ 
      names_input <- names(evt_list$time_data)
      prev_value <- setNames(vector("list", length(names_input)), names_input)
      prev_value[names_input] <- evt_list$time_data[names_input]
      prev_value["cur_evtlist"] <- list(setNames(rep(Inf,length(input_list_arm$init_event_list[[1]]$evts)), input_list_arm$init_event_list[[1]]$evts))
      dump_info <- list(
        list(
          prev_value = prev_value,
          cur_value  = c(evt_list[["time_data"]],evt_list["cur_evtlist"])
        )
      )
      
      names(dump_info) <- paste0("Analysis: ", input_list_arm$sens," ", input_list_arm$sens_name_used,
                                 "; Sim: ", input_list_arm$sim,
                                 "; Patient: ", input_list_arm$i,
                                 "; Initialize Time to Events for Patient-Arm"
      )
      
      temp_log <- c(temp_log,dump_info)
    }
  

      input_list_arm <- c(input_list_arm,evt_list$time_data,evt_list["cur_evtlist"])
      
      

      # 3 Loop per event --------------------------------------------------------
      #Main environment of reference is this one
      list_env <- list(list_env = environment())

      input_list_arm <- c(input_list_arm, list_env)
      this_patient[[arm]]$evtlist <- NULL

      input_list_arm <- c(input_list_arm,output_list)

      n_evt <- 0
      while(input_list_arm$curtime < Inf){

        # Get next event, process, repeat
        output_nxtevt <- get_next_evt(input_list_arm[["cur_evtlist"]])
        Evt <- output_nxtevt$out
        input_list_arm[['cur_evtlist']] <- output_nxtevt[["evt_list"]]

        n_evt <- n_evt +1


        if (is.null(Evt)==F){  
          
          #Evalaute event
          input_list_arm <- react_evt(Evt, arm, input_list_arm)
          
          #Get extra objects to be exported
          
          if(input_list$accum_backwards){
          extra_data <- input_list_arm[c(input_list_arm$input_out,
                                         paste0(input_list_arm$uc_lists$ongoing_inputs,"_lastupdate",recycle0 = TRUE)
                                         )]
          }else{
            extra_data <- input_list_arm[c(input_list_arm$input_out)]
          }
          extra_data <- extra_data[!sapply(extra_data,is.null)]
 
              this_patient[[arm]]$evtlist[[n_evt]] <- c(evtname = Evt$evt ,
                                                        evttime = Evt$evttime,
                                                        pat_id = i,
                                                        arm = arm,
                                                        extra_data
              )
            

        } else {input_list_arm$curtime <- Inf} #if no events, stop
        
        
        
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
      final_output$error_m <- message(e$message,
                                     message("Error in analysis: ", sens,
                                             "; simulation: ", simulation,
                                             "; patient: ", if(exists("i")){i},
                                             "; arm: ", if(exists("arm")){arm},
                                             "; event: ", if(exists("Evt")){Evt$evt},
                                             "; time: ", if(exists("Evt")){Evt$evttime})
                                     )
      return(final_output)
    }else if(input_list$continue_on_error){
      final_output <- list()  
      
      final_output$error_m <- message(e$message,
                                     message("Error in analysis: ", sens,
                                             "; simulation: ", simulation,
                                             "; patient: ", if(exists("i")){i},
                                             "; arm: ", if(exists("arm")){arm},
                                             "; event: ", if(exists("Evt")){Evt$evt},
                                             "; time: ", if(exists("Evt")){Evt$evttime})
      )
      return(final_output)
      
    } else{
      stop(message(e$message,
                   message("Error in analysis: ", sens,
                           "; simulation: ", simulation,
                           "; patient: ", if(exists("i")){i},
                           "; arm: ", if(exists("arm")){arm},
                           "; event: ", if(exists("Evt")){Evt$evt},
                           "; time: ", if(exists("Evt")){Evt$evttime})
      ))
    }
  } )

}
