#' Run the simulation engine
#'
#' @param trt_list A vector of the names of the interventions evaluated in the simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param input_list A list of all other inputs: drc,drq,psa_bool,init_event_list,evt_react_list,
#' uc_lists = list(util_ongoing_list,util_instant_list,util_cycle_list,cost_ongoing_list,cost_instant_list,cost_cycle_list),
#' input_out,ipd,trt_list,simulation,npats,n_sim
#'
#' @return A data frame with the simulation results
#' @importFrom foreach %dopar% foreach
#' @importFrom purrr map map_dbl
#' @importFrom data.table rbindlist :=
#' @importFrom doParallel stopImplicitCluster
#' 
#' @keywords internal


run_engine <- function(trt_list,
                      common_pt_inputs=NULL,
                      unique_pt_inputs=NULL,
                      input_list = NULL){
  
  # Initial set-up --------------------------
  trt_list <- trt_list
  simulation <- input_list$simulation
  n_sim <- input_list$n_sim
  npats <- input_list$npats
  psa_bool <- input_list$psa_bool

  #1 Loop per patient ----------------------------------------------------------
  patdata <- vector("list", length=npats) # empty list with npats elements
  # Outer loop, repeat for each patient
  patdata <- foreach(i = 1:npats,
                     .packages = (.packages()),
                     .export = unique(c("input_list",ls(.GlobalEnv),ls(parent.env(environment())),ls(environment()))),
                     .combine = 'c') %dopar% {

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
                         output_list <- list(curtime = 0, thslys = 0, thsqalys = 0, thscosts = 0, itemlys = 0, itemqalys = 0, itemcosts = 0, thslys_undisc = 0, thsqalys_undisc = 0, thscosts_undisc = 0, itemlys_undisc = 0, itemqalys_undisc = 0, itemcosts_undisc = 0)
                         this_patient[[trt]][["thslys"]] <- 0
                         this_patient[[trt]][["thsqalys"]]<- 0
                         this_patient[[trt]][["thscosts"]]<- 0
                         this_patient[[trt]][["thslys_undisc"]] <- 0
                         this_patient[[trt]][["thsqalys_undisc"]]<- 0
                         this_patient[[trt]][["thscosts_undisc"]]<- 0
                         
                         #Extract the inputs that are unique for each patient-intervention
                         input_list_trt <- NULL
                         input_list_trt <- c(input_list_pt,list(trt=trt))
                         if(!is.null(unique_pt_inputs)){
                           for (inp in 1:length(unique_pt_inputs)) {
                             list.unique_pt_inputs <- lapply(unique_pt_inputs[inp],function(x) eval(x, input_list_trt))
                             if (!is.null(names(list.unique_pt_inputs[[1]]))) {
                               warning("Item ", names(list.unique_pt_inputs), " is named. It is strongly advised to assign unnamed objects if they are going to be processed in the model, as they can create errors depending on how they are used within the model")
                             }
                             input_list_trt <- c(input_list_trt,list.unique_pt_inputs)

                           }
                         }
                         #Make sure there are no duplicated inputs in the model, if so, take the last one
                         duplic <- duplicated(names(input_list_trt),fromLast = T)
                         if (sum(duplic)>0) { warning("Duplicated items detected, using the last one added")  }
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

                           if (is.null(Evt)==F){ #go until there are no other events
                             
                             #Evaluate event
                             input_list_trt <- react_evt(Evt, trt, input_list_trt)
                             
                             #Get extra objects to be exported
                             extra_data <- input_list_trt[input_list_trt$input_out]
                             extra_data <- extra_data[!sapply(extra_data,is.null)]
                             
                             #Save actual event list and times
                             if (length(extra_data)==0) { #if no extra data needs to be saved, omit that step
                               this_patient[[trt]]$evtlist[[n_evt]] <-   list(evtname = Evt$evt ,
                                                                              evttime = Evt$evttime,
                                                                              pat_id = i,
                                                                              trt = trt,
                                                                              cost = input_list_trt[['itemcosts']],
                                                                              cost_undisc = input_list_trt[['itemcosts_undisc']],
                                                                              qaly = input_list_trt[['itemqalys']],
                                                                              qaly_undisc = input_list_trt[['itemqalys_undisc']],
                                                                              ly = input_list_trt[['itemlys']],
                                                                              ly_undisc = input_list_trt[['itemlys_undisc']]
                               )

                             } else{
                               this_patient[[trt]]$evtlist[[n_evt]] <- c(evtname = Evt$evt ,
                                                                         evttime = Evt$evttime,
                                                                         pat_id = i,
                                                                         trt = trt,
                                                                         cost = input_list_trt[['itemcosts']],
                                                                         cost_undisc = input_list_trt[['itemcosts_undisc']],
                                                                         qaly = input_list_trt[['itemqalys']],
                                                                         qaly_undisc = input_list_trt[['itemqalys_undisc']],
                                                                         ly = input_list_trt[['itemlys']],
                                                                         ly_undisc = input_list_trt[['itemlys_undisc']],
                                                                         extra_data
                               )
                             }
                            
                             #Accumulate total lys,costs and qalys
                             this_patient[[trt]][["thslys"]]   <- this_patient[[trt]][["thslys"]] + input_list_trt[['itemlys']]
                             this_patient[[trt]][["thsqalys"]] <- this_patient[[trt]][["thsqalys"]] + input_list_trt[['itemqalys']]
                             this_patient[[trt]][["thscosts"]] <- this_patient[[trt]][["thscosts"]]  + input_list_trt[['itemcosts']]
                             this_patient[[trt]][["thslys_undisc"]]   <- this_patient[[trt]][["thslys_undisc"]] + input_list_trt[['itemlys_undisc']]
                             this_patient[[trt]][["thsqalys_undisc"]] <- this_patient[[trt]][["thsqalys_undisc"]] + input_list_trt[['itemqalys_undisc']]
                             this_patient[[trt]][["thscosts_undisc"]] <- this_patient[[trt]][["thscosts_undisc"]]  + input_list_trt[['itemcosts_undisc']]
                             
                           } else {input_list_trt$curtime <- Inf} #if no events, stop
                           
                         }

                       }

                       return(list(this_patient))
                     }
  
  # Organize and create output -----------------------------------------------------------

final_output <- list()

#Create total outputs and user-defined costs/utilities from IPD
  vector_total_outputs <- c("total_lys_","total_qalys_","total_costs_","total_lys_undisc_","total_qalys_undisc_","total_costs_undisc_")
  vector_total_outputs_search <- c("thslys","thsqalys","thscosts","thslys_undisc","thsqalys_undisc","thscosts_undisc")
  
  vector_other_outputs <- input_list$categories_for_export
  for (trt in trt_list) {
    for (output_i in 1:length(vector_total_outputs)) {
      final_output[[paste0(vector_total_outputs[output_i],trt)]] <- sum(unlist(map(map(patdata,trt),vector_total_outputs_search[output_i])))/npats
    }
    
    for (output_i in 1:length(vector_other_outputs)) {
      final_output[[paste0(vector_other_outputs[output_i],"_",trt)]] <- sum(unlist(map_depth(map(map(patdata,trt),"evtlist"),2,vector_other_outputs[output_i])))/npats
    }
  }
  
  final_output$trt_list <- trt_list
  
#Exports IPD values
if (input_list$ipd==TRUE) {
  merged_df <- NULL
  for (trt in trt_list) {
    merged_df <- rbindlist(list(merged_df,rbindlist(unlist(map(map(patdata,trt),"evtlist"), recursive = FALSE))))
  }

  for (trt_ch in trt_list) {
    thscosts_pat <- map_dbl(map(patdata,trt_ch),"thscosts")
    thslys_pat <- map_dbl(map(patdata,trt_ch),"thslys")
    thsqalys_pat <- map_dbl(map(patdata,trt_ch),"thsqalys")
    thscosts_pat_undisc <- map_dbl(map(patdata,trt_ch),"thscosts_undisc")
    thslys_pat_undisc <- map_dbl(map(patdata,trt_ch),"thslys_undisc")
    thsqalys_pat_undisc <- map_dbl(map(patdata,trt_ch),"thsqalys_undisc")
    names(thscosts_pat) <- 1:npats
    names(thsqalys_pat) <- 1:npats
    names(thslys_pat) <- 1:npats
    names(thscosts_pat_undisc) <- 1:npats
    names(thsqalys_pat_undisc) <- 1:npats
    names(thslys_pat_undisc) <- 1:npats
    merged_df[trt==trt_ch,total_costs:= thscosts_pat[match(merged_df[trt==trt_ch,pat_id],names(thscosts_pat))]]
    merged_df[trt==trt_ch,total_costs_undisc:= thscosts_pat_undisc[match(merged_df[trt==trt_ch,pat_id],names(thscosts_pat_undisc))]]
    merged_df[trt==trt_ch,total_qalys:= thsqalys_pat[match(merged_df[trt==trt_ch,pat_id],names(thsqalys_pat))]]
    merged_df[trt==trt_ch,total_qalys_undisc:= thsqalys_pat_undisc[match(merged_df[trt==trt_ch,pat_id],names(thsqalys_pat_undisc))]]
    merged_df[trt==trt_ch,total_lys:= thslys_pat[match(merged_df[trt==trt_ch,pat_id],names(thslys_pat))]]
    merged_df[trt==trt_ch,total_lys_undisc:= thslys_pat_undisc[match(merged_df[trt==trt_ch,pat_id],names(thslys_pat_undisc))]]
  }

  final_output$merged_df <- merged_df

}
gc()

  return(final_output)

}
