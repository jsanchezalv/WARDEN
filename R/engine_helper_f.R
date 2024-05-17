# Initial event list --------------------------------------------------------------------------------------------------------------------------------------

#' Execute the initial time to events and separate the events from other inputs that are stored
#'
#' @param arm_name A character string of the name of the intervention
#' @param input_list_arm A list of simulation inputs
#'
#' @return A named vector of initial event times, and a named vector of other inputs to be stored
#'
#' @examples
#' initiate_evt(arm = "int",input_list_arm = input_list_arm)
#'
#' @keywords internal
#' @noRd


initiate_evt <- function(arm_name,input_list_arm){
  position <- which(arm_name==names(input_list_arm$init_event_list))

  time_data <- local({
    evts_v <- input_list_arm$init_event_list[[position]][["evts"]]
    
    othert_v <- input_list_arm$init_event_list[[position]][["other_inp"]]
    
    list2env(mget(c(evts_v,othert_v),ifnotfound=Inf), envir=environment()) #initialize
    
    eval(input_list_arm$init_event_list[[position]][["expr"]]) #run script
    
    evttime <- lapply(mget(evts_v,ifnotfound=Inf),unname) #get event times and make sure they are unnamed
    
    othertime <- if(!is.null(othert_v)){mget(othert_v,ifnotfound=Inf)} else{NULL}  #get other inputs times
    
    out <- list(evttime=evttime, othertime=othertime)
  },input_list_arm)

  #Event data
  cur_evtlist <- unlist(time_data$evttime)
  
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
#' get_next_evt(evt_list = input_list_arm$cur_evtlist)
#'
#' @keywords internal
#' @noRd

get_next_evt <- function(evt_list){                  # This function identifies which event is to be processed next for each patient, depending on intervention

  if (length(evt_list)>0) {
    
    # min_evt <- which.min(evt_list) #old method

    #select the position in the vector that has the minimum time, adding the priority times to make sure to select the right one,
    # 4x slower than old method
    
    min_evt <- which.min(evt_list + parent.frame()$input_list_arm$precision_times[names(evt_list)]) 
    cur_evtlist <- list(out = list(evt = names(evt_list[min_evt]), evttime = evt_list[[min_evt]]), evt_list = evt_list[-min_evt])
  } else {
    cur_evtlist <- NULL
  }


  return(cur_evtlist)
}



# Reaction to Event ---------------------------------------------------------------------------------------------------------------------------------------

#' Evaluates the reactions of the event identified by GetNextEvt
#'
#' @param thisevt A two element list containing the first list from GetNextEvt: evt and evttime
#' @param arm A character string of the name of the intervention currently being processed
#' @param input_list_arm A list of simulation inputs
#'
#' @return The updated input list with after the reaction to the event is evaluated
#'
#' @examples
#' react_evt(thisevt="evt1",arm="int",input_list_arm=input_list_arm)
#'
#' @keywords internal
#' @noRd

react_evt <- function(thisevt,arm,input_list_arm=NULL){      # This function processes the next event (as identified in the GetNextEvt function)
  # Initial set-up --------------------------
  evt <- thisevt$evt                  # Identify event type
  prevtime <- input_list_arm$curtime                 # Identify time of previous event
  curtime <- thisevt$evttime         # Identify time of next event
  
  if (curtime<prevtime) {
    stop("Time of event '", evt,"': ", round(curtime,4), " is smaller than the time of previous event '", if(!is.list(input_list_arm$evt)){input_list_arm$evt}else{NA}, "': ", round(prevtime,4), ". Arm: ", arm, ", id: ", input_list_arm$i)
  }
  
  input_list_arm[["curtime"]] <- curtime
  input_list_arm[["evt"]] <- evt
  input_list_arm[["prevtime"]] <- prevtime
  input_list_arm[["arm"]] <- arm
  
  # Create costs and utilities for event --------------------------------------------------
  evt_arm <- paste(evt,arm,sep = "_")
  
  #Reset instantaneous costs/qalys/others

  if(!is.null(input_list_arm$uc_lists$instant_inputs)){
    input_list_arm[input_list_arm$uc_lists$instant_inputs] <- 0
  }
  
  if(!is.null(input_list_arm$uc_lists$ongoing_inputs)){
    input_list_arm[paste0(input_list_arm$uc_lists$ongoing_inputs,"_lastupdate")] <- 0
  }
  

  #Evaluate the reaction to the event
  input_list_arm <- eval_reactevt(input_list_arm$evt_react_list, evt,input_list_arm)


  return(input_list_arm)

}


# Evaluate event ------------------------------------------------------------------------------------------------------------------------------------------


#' Calculates the expression which has been defined in the reaction of the event 
#'
#' @param react_list The evt_react_list from the input_list_arm object. It contains the reactions to events.
#' @param evt_name The current event being processed
#' @param input_list_arm A list of simulation inputs
#'
#' @return The modified input list with the updates after executing the corresponding reactions
#'
#' @examples
#' eval_reactevt(react_list = input_list_arm$evt_react_list,evt_name ="evt1",input_list_arm=input_list_arm)
#'
#' @keywords internal
#' @noRd

eval_reactevt <-  function(react_list,evt_name,input_list_arm=NULL){
  # Initial set-up --------------------------

  position <- which(evt_name==names(react_list))
  pos_l <- length(position)
  if (pos_l==0 | pos_l>1 ) {
    stop("Reaction to event ", evt_name, " not recognised or more than one reaction found. Make sure that only one reaction has been defined for the event")    
  }

# Evaluate reaction -------------------------------------------------------

    
    
    input_list_arm <- local({
      input_list_arm <- input_list_arm
      eval(react_list[[position]][["react"]]) #run script
      out <- input_list_arm
    },input_list_arm)
    
    
    
  
  return(input_list_arm)
  

}


# Get Inputs for Event -----------------------------------------------------------------------------------------------------------------------------------------------

#' Evaluates the cost/utility/cycle unevaluated expressions to be processed by the simulation engine
#'
#' @param x The specific cost/utility and its type (ongoing, instant...) to be used, created through add_cost/add_util
#' @param ifnull Value to be used if the input has not been defined
#' @param type Identifies what type of input is being used. Can be "cost", "util", "cycle_l" (cycle length) and "cycle_starttime" (starting time of the cycle)
#' @param evt_arm_i The event-intervention identifier to understand which specific input to use, separated by an underscore
#' @param input_list_arm_i  A list of simulation inputs
#'
#' @return A numeric vector of evaluated costs/utilities/cycle lengths/starting times for the specific event and intervention defined
#'
#' @examples
#' get_input(x = input_list_arm$uc_lists$cost_ongoing_list,ifnull=0,type="cost",evt_arm_i="evt1_int",input_list_arm_i=input_list_arm)
#'
#' @keywords internal
#' @noRd

get_input <-  function(x,ifnull=0,type,evt_arm_i =evt_arm, input_list_arm_i=input_list_arm){
  out <- NULL
  items <- x[names(x)==evt_arm_i]
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
                  eval(items[[1]][[type]],input_list_arm_i)
                } #lazy eval will give error on null, so just put 0 in that case
      )
      #otherwise do loop per item
    } else{
      for (i in 1:items_l) {
        out <-  c(out,
                  if(is.null(items[[i]][[type]])){
                    ifnull
                  } else{
                    eval(items[[i]][[type]],input_list_arm_i)
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
#' @param output_sim The output_psa data frame from the list object returned by `run_sim()`
#' @param element Variable for which mean and 95% CIs are computed (single string)
#' @param round_digit Number of digits to round outputs
#'
#' @return Mean and 95% CI from the PSA samples
#' @importFrom purrr map_dbl
#' @importFrom stats quantile
#'
#' @examples
#' interval_out(output_sim=results$output_sim[[1]],element="costs.",round_digit=3)
#'
#' @keywords internal
#' @noRd

interval_out <- function(output_sim, element,round_digit=2) {
  out <- paste0(apply(do.call(rbind, map(output_sim,element)),2,function(x) format(round(mean(x,na.rm=TRUE),round_digit), big.mark=",", scientific=FALSE)),
                " (",
                apply(do.call(rbind, map(output_sim,element)),2,function(x) format(round(quantile(x,0.025,na.rm=TRUE),round_digit), big.mark=",", scientific=FALSE)) ,
                "; ",
                apply(do.call(rbind, map(output_sim,element)),2,function(x) format(round(quantile(x,0.975,na.rm=TRUE),round_digit), big.mark=",", scientific=FALSE)) ,
                ")"
  )
  return(out)
}


# Compute and Format outputs -------------------------

#' Compute the discounting and format the outputs from the simulation
#'
#' @param patdata The list with the data from the patient-treatment iterations for a single simulation
#' @param input_list The list that contains the main inputs used for the simulation
#'
#' @return List with the outputs formatted 
#'
#' @import data.table
#' @importFrom utils tail
#' @importFrom zoo na.locf

#' 
#' @details
#' It computes the discounted and undiscounted lys/costs/qalys. 
#' 
#' It also computes the mean aggregate outcomes per arm for the simulation for numeric values of length 1, discarding NAs and Inf. 
#' Extra data defined by user of length >1 (e.g., matrix) will not be displayed as part of the final IPD data.table and 
#' will instead be reported in the `extradata_raw` list.
#' 
#' Extra data defined by user of length==1 will be integrated in the final IPD data.table. 
#' 
#' For `input_out` items that are of length ==1, using ipd = 2 in `run_sim` will take the last observation per patient.
#' If using ipd = 3, it will return the average across patients of the last observation per patient (if numeric. If not numeric, it will be discarded).
#'
#' @examples
#' compute_outputs(patdata=patdata,input_list=input_list)
#'
#' @keywords internal
#' @noRd

compute_outputs <- function(patdata,input_list) {
  arm_list <- input_list$arm_list
  simulation <- input_list$simulation
  sens <- input_list$sens
  n_sim <- input_list$n_sim
  npats <- input_list$npats
  psa_bool <- input_list$psa_bool
  
  patdata_dt <- NULL
  list_patdata <- NULL

  #Split the data as to be exported as a data.table, and the extra data the user described
  data_export_aslist <- input_list$input_out[!input_list$input_out %in% input_list$categories_for_export]
  data_export_summarized_nonumeric <- data_export_aslist
    
  for (arm_i in arm_list) {
    list_evts <- unlist(map(map(patdata,arm_i),"evtlist"), recursive = FALSE)
    list_patdata <- c(list_patdata,list_evts)
  }  

  #We exclude the extra data the user described that has a length > 1 (e.g., a matrix) from the data.table
  #as there could be matrices or other objects not suitable for data.table
  
  items_length_greater_than_one <- unlist(list_patdata,recursive=FALSE)
  items_length_greater_than_one <- items_length_greater_than_one[names(items_length_greater_than_one) %in% data_export_aslist]
  
  items_length_one_numeric <- unlist(lapply(items_length_greater_than_one, function(x) length(x)==1 & is.numeric(x)))
  items_length_one_numeric <- items_length_one_numeric[items_length_one_numeric==TRUE]
  
  items_length_greater_than_one <- unlist(lapply(items_length_greater_than_one, function(x) length(x)>1))
  items_length_greater_than_one <- items_length_greater_than_one[items_length_greater_than_one==TRUE]
  
  data_export_aslist <- unique(names(items_length_greater_than_one))
  data_export_tobesummarized <- unique(names(items_length_one_numeric))
  data_export_summarized_nonumeric <- data_export_summarized_nonumeric[!data_export_summarized_nonumeric %in% c(data_export_aslist,data_export_tobesummarized)]
  
  if (length(data_export_aslist)>0) {
    list_evts <- lapply(list_patdata,function(x) x[!names(x) %in% data_export_aslist])
  }
  
  patdata_dt <- rbindlist(list(patdata_dt,rbindlist(list_patdata,fill=TRUE)))
  
  #Extract only extra data that the user wants to export
  export_list_ipd <- lapply(list_patdata,function(x) x[data_export_aslist])
  
  
  prepared_outputs_v <- c(
    if(!is.null(input_list$uc_lists$ongoing_inputs)){as.vector(outer(input_list$uc_lists$ongoing_inputs, c("undisc"), paste, sep="_"))},
    if(!is.null(input_list$uc_lists$instant_inputs)){as.vector(outer(input_list$uc_lists$instant_inputs, c("undisc"), paste, sep="_"))},
    if(!is.null(input_list$uc_lists$cycle_inputs)){as.vector(outer(input_list$uc_lists$cycle_inputs, c("undisc"), paste, sep="_"))}
  )
  
  #Use the data.table to initialize values
  patdata_dt[,prevtime:=data.table::shift(evttime,fill=0)]
  patdata_dt[,prevtime:=ifelse(prevtime>evttime,0,prevtime)]
  cols_init <- c("lys",
                 "qalys",
                 "costs",
                 "lys_undisc",
                 "qalys_undisc",
                 "costs_undisc",
                 prepared_outputs_v
                 )
  patdata_dt[,(cols_init):=0]
  
  # Discounting of Outcomes-------------------------------------------------------------
  
  #Discount and undiscount ongoing
  
  for (cat in input_list$uc_lists$ongoing_inputs) {
    
    #Set final observation to be considered as an update
    set(patdata_dt, i = patdata_dt[,.I[.N],by=.(pat_id, arm)][['V1']], j = paste0(cat,"_lastupdate"), value = 1)
    #Updated values are kept
    patdata_dt[get(paste0(cat,"_lastupdate")) == 1, value_new := get(cat)]
    #Other values are overrwritten backwards
    patdata_dt[, value_new := zoo::na.locf(value_new,fromLast=TRUE), by=.(pat_id, arm)]
    patdata_dt[, paste0(cat) := value_new]
    patdata_dt[, value_new := NULL]
    patdata_dt[, paste0(paste0(cat,"_lastupdate")) := NULL]
    
    patdata_dt[,paste0(cat,"_","undisc") := disc_ongoing_v(lcldr=0,
                                                                        lclprvtime=prevtime,
                                                                        lclcurtime=evttime,
                                                                        lclval=get(cat))]
    
    patdata_dt[,paste0(cat) := disc_ongoing_v(lcldr=input_list$drc,
                                                                 lclprvtime=prevtime,
                                                                 lclcurtime=evttime,
                                                                 lclval=get(cat))]

    if(cat %in% input_list$uc_lists$cost_categories_ongoing){
      patdata_dt[, "costs" := costs + get(cat)]
      patdata_dt[, "costs_undisc" := costs_undisc + get(paste0(cat,"_","undisc"))]
    }
    
    if(cat %in% input_list$uc_lists$util_categories_ongoing){
      patdata_dt[, "qalys" := qalys+ get(cat)]
      patdata_dt[, "qalys_undisc" := qalys_undisc + get(paste0(cat,"_","undisc"))]
    }
    

  }
  
  
  #Discount and undiscount instant
  for (cat in input_list$uc_lists$instant_inputs) {
    patdata_dt[,paste0(cat,"_","undisc") := disc_instant_v(lcldr=0,
                                                                        lclcurtime=evttime,
                                                                        lclval=get(cat))]
    
    patdata_dt[,paste0(cat) := disc_instant_v(lcldr=input_list$drc,
                                                                 lclcurtime=evttime,
                                                                 lclval=get(cat))]
    
    if(cat %in% input_list$uc_lists$cost_categories_instant){
      patdata_dt[, "costs" := costs + get(cat)]
      patdata_dt[, "costs_undisc" := costs_undisc + get(paste0(cat,"_","undisc"))]
    }
    
    if(cat %in% input_list$uc_lists$util_categories_instant){
      patdata_dt[, "qalys" := qalys + get(cat)]
      patdata_dt[, "qalys_undisc" := qalys_undisc + get(paste0(cat,"_","undisc"))]
    }
    
  }
  
  #Discount and undiscount cycle
  for (cat in input_list$uc_lists$cycle_inputs) {
    patdata_dt[,paste0(cat,"_","undisc") := disc_cycle_v(lcldr=0,
                                                                    lclprvtime=prevtime,
                                                                    cyclelength = get(paste0(cat,"_","cycle_l")),
                                                                    lclcurtime=evttime,
                                                                    lclval= get(cat),
                                                                    starttime = get(paste0(cat,"_","cycle_starttime")))] 
    
    patdata_dt[,paste0(cat) := disc_cycle_v(lcldr=input_list$drc,
                                                             lclprvtime=prevtime,
                                                             cyclelength = get(paste0(cat,"_","cycle_l")),
                                                             lclcurtime=evttime,
                                                             lclval= get(cat),
                                                             starttime = get(paste0(cat,"_","cycle_starttime")))]
    
    if(cat %in% input_list$uc_lists$cost_categories_cycle){
      patdata_dt[, "costs" := costs + get(cat)]
      patdata_dt[, "costs_undisc" := costs_undisc + get(paste0(cat,"_","undisc"))]
    }
    
    if(cat %in% input_list$uc_lists$util_categories_cycle){
      patdata_dt[, "qalys" := qalys + get(cat)]
      patdata_dt[, "qalys_undisc" := qalys_undisc + get(paste0(cat,"_","undisc"))]
    }
    
  }
  
  
  #Discount and undiscount LYs
  patdata_dt[,"lys" := disc_ongoing_v(lcldr=input_list$drq,
                                      lclprvtime=prevtime,
                                      lclcurtime=evttime,
                                      lclval=1)]
  
  patdata_dt[,"lys_undisc" := disc_ongoing_v(lcldr=0,
                                             lclprvtime=prevtime,
                                             lclcurtime=evttime,
                                             lclval=1)]
  
  #Calculate total outcomes
  patdata_dt[,"total_costs" := sum(costs),by=.(pat_id,arm)]
  patdata_dt[,"total_qalys" := sum(qalys),by=.(pat_id,arm)]
  patdata_dt[,"total_lys" := sum(lys),by=.(pat_id,arm)]
  patdata_dt[,"total_costs_undisc" := sum(costs_undisc),by=.(pat_id,arm)]
  patdata_dt[,"total_qalys_undisc" := sum(qalys_undisc),by=.(pat_id,arm)]
  patdata_dt[,"total_lys_undisc" := sum(lys_undisc),by=.(pat_id,arm)]
  
  #Partially order the data
  data.table::setcolorder(patdata_dt,   c("evtname", "evttime", "prevtime", "pat_id", "arm",
                                          "total_lys","total_qalys","total_costs",
                                          "total_costs_undisc", "total_qalys_undisc", "total_lys_undisc",
                                          "lys","qalys","costs",
                                          "lys_undisc", "qalys_undisc",  "costs_undisc"))
  
  
  # Organize and create output -----------------------------------------------------------
  
  final_output <- list()
  
  #Create total outputs and user-defined costs/utilities from IPD
  vector_total_outputs <- c("total_lys","total_qalys","total_costs","total_lys_undisc","total_qalys_undisc","total_costs_undisc")
  vector_total_outputs_search <- c("lys","qalys","costs","lys_undisc","qalys_undisc","costs_undisc")

  #Add to final outputs the total outcomes as well as the cost/utility categories totals
  vector_other_outputs <- c(input_list$categories_for_export,prepared_outputs_v)
  for (arm_i in arm_list) {
    for (output_i in 1:length(vector_total_outputs)) {
      final_output[[vector_total_outputs[output_i]]][arm_i] <- patdata_dt[arm==arm_i,.(out=sum(get(vector_total_outputs_search[output_i]),na.rm=TRUE)),by=.(pat_id)][,mean(out,na.rm=TRUE)]
    }
    for (output_i in vector_other_outputs) {
      final_output[[output_i]][arm_i] <- patdata_dt[arm==arm_i,.(out=sum(get(output_i),na.rm=TRUE)),by=.(pat_id)][,mean(out,na.rm=TRUE)]
    }
    
    for (output_i in data_export_tobesummarized) {
        #Gets last value from patient, then average for numeric
      final_output[[output_i]][arm_i] <- patdata_dt[arm==arm_i,.(out=tail(get(output_i)*is.finite(get(output_i)),n=1,na.rm=TRUE)),by=.(pat_id)][,mean(out,na.rm=TRUE)]
    }
    
    for (output_i in data_export_summarized_nonumeric) {
      #Gets last value from patient, then average for numeric
      final_output[[output_i]][arm_i] <- patdata_dt[arm==arm_i,.(out=tail(get(output_i),n=1,na.rm=TRUE)),by=.(pat_id)][,tail(out,n=1,na.rm=TRUE)]
    }
  }
  
  final_output <- c(list(arm_list=arm_list),final_output)

  #Exports IPD values either fully IPD or aggregated for events
  if (input_list$ipd==1) {
    final_output$merged_df <- patdata_dt
    if(length(data_export_aslist)>0){
      final_output$extradata_raw <- export_list_ipd
    }
    
  } else if (input_list$ipd==2) {
    
    #Get names of columns that will be used, 
    other_cols <- c("pat_id", "arm")
    #Columns that will not be summarized (event related columns, time and total_)
    cols_to_rm <- colnames(patdata_dt)[grepl("total_",colnames(patdata_dt)) | colnames(patdata_dt) %in% c("evtname", "evttime", "prevtime")]
    patdata_dt[,number_events:=1]
    #Numeric columns only
    numeric_c <- sapply(patdata_dt,is.numeric)
    #Columns to sum must be in the right list and also be numeric
    cols_to_sum <- colnames(patdata_dt)[!colnames(patdata_dt) %in% c(other_cols,cols_to_rm,data_export_tobesummarized) & numeric_c]
    #Other columns are left as is (takes last value of those)
    cols_to_leave_as_is <- colnames(patdata_dt)[!colnames(patdata_dt) %in% c(other_cols,cols_to_rm) & !colnames(patdata_dt) %in% c(cols_to_sum)]
    
    #Summarize the data as relevant
    patdata_temp <- patdata_dt[, lapply(.SD, sum, na.rm=TRUE), by=other_cols, .SDcols=cols_to_sum] #sum numeric variables in list
    
    if (length(cols_to_leave_as_is!=0)) {
      patdata_temp2 <- patdata_dt[, tail(.SD, 1, na.rm=TRUE), by=other_cols, .SDcols=cols_to_leave_as_is] #get last observation if other
      final_output$merged_df <- merge(patdata_temp,patdata_temp2)
    } else{
      final_output$merged_df <- patdata_temp
    }
    

    if (sens==1 & simulation==1) {
      message("Patient-arm data aggregated across events by selecting the last value for input_out items.")
      }
    
    if(length(data_export_aslist)>0){
      final_output$extradata_raw <- export_list_ipd
    }
  } else{
    if (sens==1 & simulation==1) {
      message("Data aggregated across events and patients by selecting the last value for input_out numeric items and then averaging across patients. Only last value of non-numeric items is displayed .Items with length > 1 have been discarded. ")
    }
  }
  
  return(final_output)
}
