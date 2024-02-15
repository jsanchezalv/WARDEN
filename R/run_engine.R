#' Run the simulation engine without using a parallel engine, facilitating debugging
#'
#' @param trt_list A vector of the names of the interventions evaluated in the simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param input_list A list of all other inputs: drc, drq, psa_bool, init_event_list, evt_react_list,
#' uc_lists = list(util_ongoing_list,util_instant_list,util_cycle_list,cost_ongoing_list,cost_instant_list,cost_cycle_list),
#' input_out,ipd,trt_list,simulation,npats,n_sim
#'
#' @return A data frame with the simulation results
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom data.table rbindlist
#'
#' @keywords internal


run_engine <- function(trt_list,
                            common_pt_inputs=NULL,
                            unique_pt_inputs=NULL,
                            input_list = NULL){
  # Initial set-up --------------------------
  trt_list <- trt_list
  simulation <- input_list$simulation
  sens <- input_list$sens
  n_sim <- input_list$n_sim
  npats <- input_list$npats
  psa_bool <- input_list$psa_bool


  #1 Loop per patient ----------------------------------------------------------
  patdata <- vector("list", length=npats) # empty list with npats elements

  for (i in 1:npats) {

    #Create empty pat data for each trt
    this_patient <- list()
    input_list_pt <- c(input_list,list(i=i))
    
    #Extract the inputs that are common for each patient across interventions
    if(!is.null(common_pt_inputs)){
      for (inp in 1:length(common_pt_inputs)) {
        list.common_pt_inputs <- lapply(common_pt_inputs[inp],function(x) eval(x, input_list_pt))
        if (!is.null(names(list.common_pt_inputs[[1]]))) {
          warning("Item ", names(list.common_pt_inputs), " is named. It is strongly advised to assign unnamed objects if they are going to be processed in the model, as they can create errors depending on how they are used within the model")
        }
        input_list_pt <- c(input_list_pt,list.common_pt_inputs)
      }
    }

    #Make sure there are no duplicated inputs in the model, if so, take the last one
    duplic <- duplicated(names(input_list_pt),fromLast = T)
    if (sum(duplic)>0) { warning("Duplicated items detected, using the last one added")  }
    input_list_pt <- input_list_pt[!duplic]

    #2 Loop per treatment ------------------------------------------------------

    for (trt in trt_list) {
      # Initialize values to prevent errors
      output_list <- list(curtime = 0)
      
      #Extract the inputs that are unique for each patient-intervention
      input_list_trt <- NULL
      input_list_trt <- c(input_list_pt,list(trt=trt))
      if(!is.null(unique_pt_inputs)){
        for (inp in 1:length(unique_pt_inputs)) {
          list.unique_pt_inputs <- lapply(unique_pt_inputs[inp],function(x) eval(x, input_list_trt))
          if ((!is.null(names(list.unique_pt_inputs[[1]]))) & i==1 & simulation==1 & sens==1) {
            warning("Item ", names(list.unique_pt_inputs), " is named. It is strongly advised to assign unnamed objects if they are going to be processed in the model, as they can create errors depending on how they are used within the model")
          }
          input_list_trt <- c(input_list_trt,list.unique_pt_inputs)

        }
      }

      #Make sure there are no duplicated inputs in the model, if so, take the last one
      duplic <- duplicated(names(input_list_trt),fromLast = T)
      if (sum(duplic)>0 & i==1 & simulation==1 & sens==1) { warning("Duplicated items detected, using the last one added")  }
      input_list_trt <- input_list_trt[!duplic]

      # Generate event list
      #if noeventlist, then just make start at 0
      if (is.null(input_list_trt$init_event_list)) {
        evt_list <- list(cur_evtlist = setNames(0,"start"), time_data = NULL)
      } else{
        evt_list <- do.call("initiate_evt",list(trt,input_list_trt))
      }


      input_list_trt <- c(input_list_trt,evt_list$time_data,evt_list["cur_evtlist"])

      # 3 Loop per event --------------------------------------------------------
      #Main environment of reference is this one
      list_env <- list(list_env = environment())

      input_list_trt <- c(input_list_trt, list_env)
      this_patient[[trt]]$evtlist <- NULL

      input_list_trt <- c(input_list_trt,output_list)

      n_evt <- 0
      while(input_list_trt$curtime < Inf){

        # Get next event, process, repeat
        output_nxtevt <- get_next_evt(input_list_trt[["cur_evtlist"]])
        Evt <- output_nxtevt$out
        input_list_trt[['cur_evtlist']] <- output_nxtevt[["evt_list"]]

        n_evt <- n_evt +1


        if (is.null(Evt)==F){  
          
          #Evalaute event
          input_list_trt <- react_evt(Evt, trt, input_list_trt)
          
          #Get extra objects to be exported
          extra_data <- input_list_trt[input_list_trt$input_out]
          extra_data <- extra_data[!sapply(extra_data,is.null)]
 
              this_patient[[trt]]$evtlist[[n_evt]] <- c(evtname = Evt$evt ,
                                                        evttime = Evt$evttime,
                                                        pat_id = i,
                                                        trt = trt,
                                                        extra_data
              )
            

        } else {input_list_trt$curtime <- Inf} #if no events, stop
        
        
        
      }

    }

    patdata[[i]] <- this_patient
  }

# Compute outputs ---------------------------------------------------------


  #Compute the outputs and format the data
  final_output <- compute_outputs(patdata, input_list)
  

  return(final_output)

}
