# Summary for deterministic/last created output for specific treatment ---------------------------

#' Deterministic results for a specific treatment
#'
#' @param out The final_output data frame from the list object returned by `run_sim()`
#' @param arm The reference treatment for calculation of incremental outcomes
#'
#' @return A dataframe with absolute costs, LYs, QALYs, and ICER and ICUR for each intervention
#' @export
#'
#' @examples
#' \dontrun{
#' summary_results_det(results$output_sim[[1]][[1]],arm="int")
#' }

summary_results_det <- function(out = results$output_sim[[1]][[1]], arm=NULL){
  arm <- ifelse(is.null(arm),out$arm_list[1],arm)

  other_arm_list <- out$arm_list[out$arm_list!=arm] #For any other treatment that is not the reference one, add the reference arm costs/lys/qalys

  remove_outputs_list <- c("total_lys_int", "total_qalys_int", "total_costs_int", "total_lys_undisc_int", 
                             "total_qalys_undisc_int", "total_costs_undisc_int", "total_lys_noint", "total_qalys_noint", 
                             "total_costs_noint", "total_lys_undisc_noint", "total_qalys_undisc_noint", 
                             "total_costs_undisc_noint", 
                             "arm_list", "merged_df")

  other_outputs <- names(out)[!names(out)%in% remove_outputs_list]
  other_outputs <- unique(sub("_[^_]+$", "", other_outputs)) #remove treatment indicator
  
  for (other_arm in other_arm_list) { #add reference arm costs/lys/qalys to compare with the other arms
    out[[paste0("dlys_",other_arm)]] <-  out[[paste0("total_lys_",arm)]] - out[[paste0("total_lys_",other_arm)]]
    out[[paste0("dqalys_",other_arm)]] <-   out[[paste0("total_qalys_",arm)]] - out[[paste0("total_qalys_",other_arm)]]
    out[[paste0("dcosts_",other_arm)]] <-   out[[paste0("total_costs_",arm)]] - out[[paste0("total_costs_",other_arm)]]
    out[[paste0("dlys_undisc_",other_arm)]] <-  out[[paste0("total_lys_undisc_",arm)]] - out[[paste0("total_lys_undisc_",other_arm)]]
    out[[paste0("dqalys_undisc_",other_arm)]] <-   out[[paste0("total_qalys_undisc_",arm)]] - out[[paste0("total_qalys_undisc_",other_arm)]]
    out[[paste0("dcosts_undisc_",other_arm)]] <-   out[[paste0("total_costs_undisc_",arm)]] - out[[paste0("total_costs_undisc_",other_arm)]]

    out[[paste0("ICER_",other_arm)]] <-   out[[paste0("dcosts_",other_arm)]] / out[[paste0("dlys_",other_arm)]]
    out[[paste0("ICUR_",other_arm)]] <-   out[[paste0("dcosts_",other_arm)]] / out[[paste0("dqalys_",other_arm)]]
    out[[paste0("ICER_undisc_",other_arm)]] <-   out[[paste0("dcosts_undisc_",other_arm)]] / out[[paste0("dlys_undisc_",other_arm)]]
    out[[paste0("ICUR_undisc_",other_arm)]] <-   out[[paste0("dcosts_undisc_",other_arm)]] / out[[paste0("dqalys_undisc_",other_arm)]]
    
    for (output_i in other_outputs) {
      out[[paste0("d",output_i,other_arm)]] <-  out[[paste0(output_i,"_",arm)]] - out[[paste0(output_i,"_",other_arm)]]
      
    }
    

  }

  out[[paste0("dlys_",arm)]] <-   0
  out[[paste0("dqalys_",arm)]] <-   0
  out[[paste0("dcosts_",arm)]] <-   0
  out[[paste0("dlys_undisc_",arm)]] <-   0
  out[[paste0("dqalys_undisc_",arm)]] <-   0
  out[[paste0("dcosts_undisc_",arm)]] <-   0  
  for (output_i in other_outputs) {
    out[[paste0("d",output_i,arm)]] <-  0
  }

  out[[paste0("ICER_",arm)]] <-   NA
  out[[paste0("ICUR_",arm)]] <-   NA
  out[[paste0("ICER_undisc_",arm)]] <-   NA
  out[[paste0("ICUR_undisc_",arm)]] <-   NA


data <- data.frame()
  for (arm in out$arm_list) {
    temp <- data.frame(
      arm = arm,
      costs = out[[paste0("total_costs_",arm)]],
      lys = out[[paste0("total_lys_",arm)]],
      qalys = out[[paste0("total_qalys_",arm)]],
      ICER = out[[paste0("ICER_",arm)]],
      ICUR = out[[paste0("ICUR_",arm)]],
      costs_undisc = out[[paste0("total_costs_undisc_",arm)]],
      lys_undisc = out[[paste0("total_lys_undisc_",arm)]],
      qalys_undisc = out[[paste0("total_qalys_undisc_",arm)]],
      ICER_undisc = out[[paste0("ICER_undisc_",arm)]],
      ICUR_undisc = out[[paste0("ICUR_undisc_",arm)]]
    )
    for (output_i in other_outputs) {
      temp[[output_i]] <-  out[[paste0(output_i,"_",arm)]]
    }
    
    data <- rbind(data,temp)
  }

  names <-  data[,1]

  # Transpose everything other than the first column
  data <- as.data.frame(as.matrix(t(data[,-1])))

  # Assign first column as the column names of the transposed dataframe
  colnames(data) <- names

  data <-round(data,2)

  return(data)

}

# Summary for PSA output for specific treatment --------------------------------------------------

#' Summary of PSA outputs for a treatment
#'
#' @param out The output_sim data frame from the list object returned by `run_sim()`
#' @param arm The reference treatment for calculation of incremental outcomes
#'
#' @return A data frame with mean and 95% CI of absolute costs, LYs, QALYs, ICER and ICUR for each intervention from the PSA samples
#' @export
#'
#' @examples
#' \dontrun{
#' summary_results_sim(results$output_sim[[1]], arm="int")
#' }

summary_results_sim <- function(out = results$output_sim[[1]], arm=NULL){

  arm <- ifelse(is.null(arm),out[[1]]$arm_list[1],arm)

  other_arm_list <- out[[1]]$arm_list[out[[1]]$arm_list!=arm] #For any other treatment that is not the reference one, add the reference arm costs/lys/qalys
  
  remove_outputs_list <- c("total_lys_int", "total_qalys_int", "total_costs_int", "total_lys_undisc_int", 
                           "total_qalys_undisc_int", "total_costs_undisc_int", "total_lys_noint", "total_qalys_noint", 
                           "total_costs_noint", "total_lys_undisc_noint", "total_qalys_undisc_noint", 
                           "total_costs_undisc_noint", 
                           "arm_list", "merged_df")
  
  other_outputs <- names(out[[1]])[!names(out[[1]])%in% remove_outputs_list]
  other_outputs <- unique(sub("_[^_]+$", "", other_outputs)) #remove treatment indicator
  other_outputs <- other_outputs[!other_outputs=="extradata"]
  
  for (sim in 1:length(out)) {

    for (other_arm in other_arm_list) { #add reference arm costs/lys/qalys to compare with the other arms
      out[[sim]][[paste0("dlys_",other_arm)]] <-  out[[sim]][[paste0("total_lys_",arm)]] - out[[sim]][[paste0("total_lys_",other_arm)]]
      out[[sim]][[paste0("dqalys_",other_arm)]] <-   out[[sim]][[paste0("total_qalys_",arm)]] - out[[sim]][[paste0("total_qalys_",other_arm)]]
      out[[sim]][[paste0("dcosts_",other_arm)]] <-   out[[sim]][[paste0("total_costs_",arm)]] - out[[sim]][[paste0("total_costs_",other_arm)]]
      out[[sim]][[paste0("dlys_undisc_",other_arm)]] <-  out[[sim]][[paste0("total_lys_undisc_",arm)]] - out[[sim]][[paste0("total_lys_undisc_",other_arm)]]
      out[[sim]][[paste0("dqalys_undisc_",other_arm)]] <-   out[[sim]][[paste0("total_qalys_undisc_",arm)]] - out[[sim]][[paste0("total_qalys_undisc_",other_arm)]]
      out[[sim]][[paste0("dcosts_undisc_",other_arm)]] <-   out[[sim]][[paste0("total_costs_undisc_",arm)]] - out[[sim]][[paste0("total_costs_undisc_",other_arm)]]
      
      out[[sim]][[paste0("ICER_",other_arm)]] <-   out[[sim]][[paste0("dcosts_",other_arm)]] / out[[sim]][[paste0("dlys_",other_arm)]]
      out[[sim]][[paste0("ICUR_",other_arm)]] <-   out[[sim]][[paste0("dcosts_",other_arm)]] / out[[sim]][[paste0("dqalys_",other_arm)]]
      out[[sim]][[paste0("ICER_undisc_",other_arm)]] <-   out[[sim]][[paste0("dcosts_undisc_",other_arm)]] / out[[sim]][[paste0("dlys_undisc_",other_arm)]]
      out[[sim]][[paste0("ICUR_undisc_",other_arm)]] <-   out[[sim]][[paste0("dcosts_undisc_",other_arm)]] / out[[sim]][[paste0("dqalys_undisc_",other_arm)]]
      
      for (output_i in other_outputs) {
        out[[sim]][[paste0("d",output_i,other_arm)]] <-  out[[sim]][[paste0(output_i,"_",arm)]] - out[[sim]][[paste0(output_i,"_",other_arm)]]
        
      }
      
    }

    out[[sim]][[paste0("dlys_",arm)]] <-   0
    out[[sim]][[paste0("dqalys_",arm)]] <-   0
    out[[sim]][[paste0("dcosts_",arm)]] <-   0
    out[[sim]][[paste0("dlys_undisc_",arm)]] <-   0
    out[[sim]][[paste0("dqalys_undisc_",arm)]] <-   0
    out[[sim]][[paste0("dcosts_undisc_",arm)]] <-   0  
    for (output_i in other_outputs) {
      out[[sim]][[paste0("d",output_i,arm)]] <-  0
    }
    
    out[[sim]][[paste0("ICER_",arm)]] <-   NA
    out[[sim]][[paste0("ICUR_",arm)]] <-   NA
    out[[sim]][[paste0("ICER_undisc_",arm)]] <-   NA
    out[[sim]][[paste0("ICUR_undisc_",arm)]] <-   NA

  }

  data <- data.frame()
  for (arm in out[[1]]$arm_list) {
    temp <- data.frame(
      arm = arm,
      costs = interval_out(out,"total_costs_",arm,0),
      lys = interval_out(out,"total_lys_",arm,2),
      qalys =  interval_out(out,"total_qalys_",arm,2),
      ICER =  interval_out(out,"ICER_",arm,0),
      ICUR =  interval_out(out,"ICUR_",arm,0),
      costs_undisc = interval_out(out,"total_costs_undisc_",arm,0),
      lys_undisc = interval_out(out,"total_lys_undisc_",arm,2),
      qalys_undisc =  interval_out(out,"total_qalys_undisc_",arm,2),
      ICER_undisc =  interval_out(out,"ICER_undisc_",arm,0),
      ICUR_undisc =  interval_out(out,"ICUR_undisc_",arm,0)
    )
    
    for (output_i in other_outputs) {
      temp[[output_i]] <-  interval_out(out,paste0(output_i,"_"),arm,2) 
    }

    data <- rbind(data,temp)
  }

  names <-  data[,1]

  # Transpose everything other than the first column
  data <- as.data.frame(as.matrix(t(data[,-1])))

  # Assign first column as the column names of the transposed dataframe
  colnames(data) <- names

  return(data)

}


# Extract all specific PSA result -------------------------------------------------------------------------------------------------------------------------

#' Extract PSA results from a treatment
#'
#' @param x The output_sim data frame from the list object returned by `run_sim()`
#' @param element Variable for which PSA results are being extracted (single string)
#' @param arm Intervention for which PSA results are being extracted (single string)
#'
#' @return A dataframe with PSA results from the specified intervention
#' @export
#'
#' @examples
#' \dontrun{
#' extract_psa_result(results$output_sim[[1]],"costs","int")
#' }

extract_psa_result <- function(x, element,arm) {
  out <- purrr::map_dbl(x,paste0(paste0(element,"_"),arm))

  out <- data.frame(element = element, arm = arm , simulation = 1:length(out),value=out)
  return(out)
}





# CEAC -------------------------------------------------------------------------------------------------------------------------

#' Calculate the cost-effectiveness acceptability curve (CEAC) for a DES model with a PSA result
#'
#' @param wtp Vector of length >=1 with the willingness to pay
#' @param results The list object returned by `run_sim()`
#' @param interventions A character vector with the names of the interventions to be used for the analysis
#' @param sensitivity_used Integer signaling which sensitivity analysis to use
#'
#' @return A data frame with the CEAC results
#' @export
#'
#' @examples
#' \dontrun{
#' ceac_des(seq(from=10000,to=500000,by=10000),results)
#' }

ceac_des <- function(wtp, results, interventions = NULL, sensitivity_used = 1) {

  if (is.null(interventions)) {
    interventions <- results$output_sim[[sensitivity_used]][[1]]$arm_list
  }

  nmb <- data.frame()
  for (comparator in interventions) {

     nmb_i <- t(as.matrix(sapply(wtp, function(wtp_i) wtp_i * extract_psa_result(results$output_sim[[sensitivity_used]],"total_qalys",comparator)$value - extract_psa_result(results$output_sim[[sensitivity_used]],"total_costs",comparator)$value)))

     if (nrow(nmb_i)>1) {
       nmb_i <- t(nmb_i)
     }
     
     nmb_i <- data.frame(nmb_i,stringsAsFactors = FALSE)
     
     names(nmb_i) <- format(wtp, scientific=F)
     nmb_i$iteration <-  1:nrow(nmb_i)
     nmb_i <- nmb_i %>% gather(key="wtp",value="nmb",-iteration)
     nmb_i$comparator <- comparator
     nmb_i$wtp <- as.numeric(nmb_i$wtp)

    nmb <- rbind(nmb,nmb_i)
  }

  nmb <- nmb %>%
    group_by(wtp,iteration) %>%
    mutate(best_nmb= ifelse(max(nmb)>0,comparator[nmb==max(nmb)],NA))

  ceac <- nmb %>%
    group_by(wtp,comparator) %>%
    summarise(prob_best= sum(best_nmb==comparator)/n()) %>%
    mutate(prob_best = ifelse(is.na(prob_best),0,prob_best))


  return(ceac)
}


# EVPI -------------------------------------------------------------------------------------------------------------------------

#' Calculate the Expected Value of Perfect Information (EVPI) for a DES model with a PSA result
#'
#' @param wtp Vector of length >=1 with the willingness to pay
#' @param results The list object returned by `run_sim()`
#' @param interventions A character vector with the names of the interventions to be used for the analysis
#' @param sensitivity_used Integer signaling which sensitivity analysis to use
#'
#' @return A data frame with the EVPI results
#' @export
#'
#' @examples
#' \dontrun{
#' evpi_des(seq(from=10000,to=500000,by=10000),results)
#' }

evpi_des <- function(wtp, results, interventions = NULL, sensitivity_used = 1) {

  if (is.null(interventions)) {
    interventions <- results$output_sim[[sensitivity_used]][[1]]$arm_list
  }


  nmb <- data.frame()
  for (comparator in interventions) {

    nmb_i <- t(as.matrix(sapply(wtp, function(wtp_i) wtp_i * extract_psa_result(results$output_sim[[sensitivity_used]],"total_qalys",comparator)$value - extract_psa_result(results$output_sim[[sensitivity_used]],"total_costs",comparator)$value)))

    if (nrow(nmb_i)>1) {
      nmb_i <- t(nmb_i)
    }
    nmb_i <- data.frame(nmb_i, stringsAsFactors = FALSE)
    
    names(nmb_i) <- format(wtp, scientific=F)
    nmb_i$iteration <-  1:nrow(nmb_i)
    nmb_i <- nmb_i %>% gather(key="wtp",value="nmb",-iteration)
    nmb_i$comparator <- comparator
    nmb_i$wtp <- as.numeric(nmb_i$wtp)

    nmb <- rbind(nmb,nmb_i)
  }

  nmb <-nmb %>%
    group_by(wtp,comparator) %>%
    mutate(mean_nmb=mean(nmb)) %>%
    group_by(wtp,iteration) %>%
    mutate(max_nmb=max(nmb)) %>%
    group_by(wtp) %>%
    mutate(max_mean_nmb = max(mean_nmb),
           mean_max_nmb = mean(max_nmb)) %>%
    ungroup() %>%
    mutate(evpi = mean_max_nmb - max_mean_nmb) %>%
    select(wtp,evpi)


  return(nmb)
}
