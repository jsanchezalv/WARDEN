
# Global variables for CRAN check -----------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      c('iteration',
        'best_nmb',
        'n',
        'prob_best',
        'results',
        'iteration',
        'mean_nmb',
        'max_nmb',
        'mean_max_nmb',
        'max_mean_nmb',
        'evpi')
    )) 
}


# Summary for deterministic/last created output for specific treatment ---------------------------

#' Deterministic results for a specific treatment
#'
#' @param out The final_output data frame from the list object returned by `run_sim()`
#' @param arm The reference treatment for calculation of incremental outcomes
#' @param wtp Willingness to pay to have INMB
#'
#' @return A dataframe with absolute costs, LYs, QALYs, and ICER and ICUR for each intervention
#' @export
#'
#' @examples
#' 
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' 
#' summary_results_det(res[[1]][[1]],arm="int")
#' 

summary_results_det <- function(out = results[[1]][[1]], arm = NULL, wtp = 50000){
  arm_ref <- ifelse(is.null(arm),1,match(arm,out$arm_list))

  remove_outputs_list <- c("total_lys", "total_qalys", "total_costs", "total_lys_undisc", 
                             "total_qalys_undisc", "total_costs_undisc",
                             "arm_list", "merged_df","sensitivity_name")
  
  outputs_names <- c("dlys","dqalys","dcosts","dlys_undisc","dqalys_undisc","dcosts_undisc")

  other_outputs <- names(out)[!names(out)%in% remove_outputs_list]
  other_outputs <- other_outputs[!other_outputs=="extradata_raw"]
  other_outputs <- other_outputs[unlist(sapply(out[other_outputs],is.numeric))]
  
  for (arm_i in 1:length(out$arm_list)) { #add reference arm costs/lys/qalys to compare with the other arms
    
    for (output in 1:length(outputs_names)) { #standard outputs
      out[[outputs_names[output]]][arm_i] <- out[[remove_outputs_list[output]]][arm_ref] - out[[remove_outputs_list[output]]][arm_i]
    }
    
    for (output_i in other_outputs) { #other outputs
      out[[paste0("d",output_i)]][arm_i] <-  out[[output_i]][arm_ref] - out[[output_i]][arm_i]
    }
    
    if (arm_i!=arm_ref) {
      out[["ICER"]][arm_i]         <-   out[["dcosts"]][arm_i]  / out[["dlys"]][arm_i] 
      out[["ICUR"]][arm_i]         <-   out[["dcosts"]][arm_i]  / out[["dqalys"]][arm_i] 
      out[["ICER_undisc"]][arm_i]  <-   out[["dcosts_undisc"]][arm_i]  / out[["dlys_undisc"]][arm_i] 
      out[["ICUR_undisc"]][arm_i]  <-   out[["dcosts_undisc"]][arm_i]  / out[["dqalys_undisc"]][arm_i] 
      out[["INMB"]][arm_i]         <-   wtp *  out[["dqalys"]][arm_i] - out[["dcosts"]][arm_i] 
      out[["INMB_undisc"]][arm_i]  <-  wtp *  out[["dqalys_undisc"]][arm_i] - out[["dcosts_undisc"]][arm_i] 
    } else{
      out[["ICER"]][arm_i]         <- NA
      out[["ICUR"]][arm_i]         <- NA
      out[["ICER_undisc"]][arm_i]  <- NA
      out[["ICUR_undisc"]][arm_i]  <- NA
      out[["INMB"]][arm_i]         <- NA
      out[["INMB_undisc"]][arm_i]  <- NA
    }
  }

  data <- data.frame(
      arm = out$arm_list,
      costs = out[["total_costs"]],
      dcosts = out[["dcosts"]],
      lys = out[["total_lys"]],
      dlys = out[["dlys"]],
      qalys = out[["total_qalys"]],
      dqalys = out[["dqalys"]],
      ICER = out[["ICER"]],
      ICUR = out[["ICUR"]],
      INMB = out[["INMB"]],
      costs_undisc = out[["total_costs_undisc"]],
      dcosts_undisc = out[["dcosts_undisc"]],
      lys_undisc = out[["total_lys_undisc"]],
      dlys_undisc = out[["dlys_undisc"]],
      qalys_undisc = out[["total_qalys_undisc"]],
      dqalys_undisc = out[["dqalys_undisc"]],
      ICER_undisc = out[["ICER_undisc"]],
      ICUR_undisc = out[["ICUR_undisc"]],
      INMB_undisc = out[["INMB_undisc"]]
    )
    for (output_i in other_outputs) {
      data[[output_i]] <-  out[[output_i]]
      data[[paste0("d",output_i)]] <-  out[[paste0("d",output_i)]]
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
#' @param wtp Willingness to pay to have INMB
#'
#' @return A data frame with mean and 95% CI of absolute costs, LYs, QALYs, ICER and ICUR for each intervention from the PSA samples
#' @export
#'
#' @examples
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' 
#' summary_results_sim(res[[1]],arm="int")

summary_results_sim <- function(out = results[[1]], arm=NULL, wtp = 50000){

  arm_ref <- ifelse(is.null(arm),1,match(arm,out[[1]]$arm_list))

  remove_outputs_list <- c("total_lys", "total_qalys", "total_costs", "total_lys_undisc", 
                           "total_qalys_undisc", "total_costs_undisc",
                           "arm_list", "merged_df","sensitivity_name")
  
  outputs_names <- c("dlys","dqalys","dcosts","dlys_undisc","dqalys_undisc","dcosts_undisc")
  
  other_outputs <- names(out[[1]])[!names(out[[1]])%in% remove_outputs_list]
  other_outputs <- other_outputs[!other_outputs=="extradata_raw"]
  other_outputs <- other_outputs[unlist(sapply(out[[1]][other_outputs],is.numeric))]
  
  for (sim in 1:length(out)) {
    if(is.null(out[[sim]])){next}

    for (arm_i in 1:length(out[[1]]$arm_list)) { #add reference arm costs/lys/qalys to compare with the other arms
      
      for (output in 1:length(outputs_names)) { #standard outputs
        out[[sim]][[outputs_names[output]]][arm_i] <- out[[sim]][[remove_outputs_list[output]]][arm_ref] - out[[sim]][[remove_outputs_list[output]]][arm_i]
      }
      
      for (output_i in other_outputs) { #other outputs
        out[[sim]][[paste0("d",output_i)]][arm_i] <-  out[[sim]][[output_i]][arm_ref] - out[[sim]][[output_i]][arm_i]
      }
      
      if (arm_i!=arm_ref) {
        out[[sim]][["ICER"]][arm_i]         <-   out[[sim]][["dcosts"]][arm_i]  / out[[sim]][["dlys"]][arm_i] 
        out[[sim]][["ICUR"]][arm_i]         <-   out[[sim]][["dcosts"]][arm_i]  / out[[sim]][["dqalys"]][arm_i] 
        out[[sim]][["ICER_undisc"]][arm_i]  <-   out[[sim]][["dcosts_undisc"]][arm_i]  / out[[sim]][["dlys_undisc"]][arm_i] 
        out[[sim]][["ICUR_undisc"]][arm_i]  <-   out[[sim]][["dcosts_undisc"]][arm_i]  / out[[sim]][["dqalys_undisc"]][arm_i] 
        out[[sim]][["INMB"]][arm_i]         <-   wtp * out[[sim]][["dqalys"]][arm_i] - out[[sim]][["dcosts"]][arm_i] 
        out[[sim]][["INMB_undisc"]][arm_i]  <-   wtp * out[[sim]][["dqalys_undisc"]][arm_i] - out[[sim]][["dcosts_undisc"]][arm_i] 
      } else{
        out[[sim]][["ICER"]][arm_i]         <- NA
        out[[sim]][["ICUR"]][arm_i]         <- NA
        out[[sim]][["ICER_undisc"]][arm_i]  <- NA
        out[[sim]][["ICUR_undisc"]][arm_i]  <- NA
        out[[sim]][["INMB"]][arm_i]         <- NA
        out[[sim]][["INMB_undisc"]][arm_i]  <- NA
      }
    }

  }

  data <- data.frame(
    arm           = out[[1]]$arm_list,
    costs         = interval_out(out,"total_costs",0),
    dcosts        = interval_out(out,"dcosts",0),
    lys           = interval_out(out,"total_lys",2),
    dlys          = interval_out(out,"dlys",3),
    qalys         = interval_out(out,"total_qalys",2),
    dqalys        = interval_out(out,"dqalys",3),
    ICER          = interval_out(out,"ICER",0),
    ICUR          = interval_out(out,"ICUR",0),
    INMB          = interval_out(out,"INMB",0),
    costs_undisc  = interval_out(out,"total_costs_undisc",0),
    dcosts_undisc = interval_out(out,"dcosts_undisc",0),
    lys_undisc    = interval_out(out,"total_lys_undisc",2),
    dlys_undisc   = interval_out(out,"dlys_undisc",3),
    qalys_undisc  = interval_out(out,"total_qalys_undisc",2),
    dqalys_undisc = interval_out(out,"dqalys_undisc",3),
    ICER_undisc   = interval_out(out,"ICER_undisc",0),
    ICUR_undisc   = interval_out(out,"ICUR_undisc",0),
    INMB_undisc   = interval_out(out,"INMB_undisc",0)
  )
  for (output_i in other_outputs) {
    data[[output_i]] <- interval_out(out,output_i,2)
    data[[paste0("d",output_i)]] <-  interval_out(out,paste0("d",output_i),3)
  }

  names <-  data[,1]

  # Transpose everything other than the first column
  data <- as.data.frame(as.matrix(t(data[,-1])))

  # Assign first column as the column names of the transposed dataframe
  colnames(data) <- names

  return(data)

}


# Summary for sensitivity output for specific treatment --------------------------------------------------

#' Summary of sensitivity outputs for a treatment
#'
#' @param out The list object returned by `run_sim()`
#' @param arm The reference treatment for calculation of incremental outcomes
#' @param wtp Willingness to pay to have INMB
#'
#' @return A data frame with each sensitivity output per arm
#' @export
#'
#' @examples
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' 
#' summary_results_sens(res,arm="int")

summary_results_sens <- function(out = results, arm=NULL, wtp = 50000){
  
  arm_ref <- ifelse(is.null(arm),1,match(arm,out[[1]][[1]]$arm_list))
  
  remove_outputs_list <- c("total_lys", "total_qalys", "total_costs", "total_lys_undisc", 
                           "total_qalys_undisc", "total_costs_undisc",
                           "arm_list", "merged_df","sensitivity_name")
  
  outputs_names <- c("dlys","dqalys","dcosts","dlys_undisc","dqalys_undisc","dcosts_undisc")
  
  other_outputs <- names(out[[1]][[1]])[!names(out[[1]][[1]])%in% remove_outputs_list]
  other_outputs <- other_outputs[!other_outputs=="extradata_raw"]
  other_outputs <- other_outputs[unlist(sapply(out[[1]][[1]][other_outputs],is.numeric))]
  
  
  data_final <- data.frame()
  
  for (sens in 1:length(out)) {
    if(is.null(out[[sens]])){next}
    data <- data.frame()
    for (sim in 1:length(out[[1]])) {
      if(is.null(out[[sens]][[sim]])){next}
      
      
      for (arm_i in 1:length(out[[1]][[1]]$arm_list)) { #add reference arm costs/lys/qalys to compare with the other arms
        
        for (output in 1:length(outputs_names)) { #standard outputs
          out[[sens]][[sim]][[outputs_names[output]]][arm_i] <- out[[sens]][[sim]][[remove_outputs_list[output]]][arm_ref] - out[[sens]][[sim]][[remove_outputs_list[output]]][arm_i]
        }
        
        for (output_i in other_outputs) { #other outputs
          out[[sens]][[sim]][[paste0("d",output_i)]][arm_i] <-  out[[sens]][[sim]][[output_i]][arm_ref] - out[[sens]][[sim]][[output_i]][arm_i]
        }
        
        if (arm_i!=arm_ref) {
          out[[sens]][[sim]][["ICER"]][arm_i]         <-   out[[sens]][[sim]][["dcosts"]][arm_i]             / out[[sens]][[sim]][["dlys"]][arm_i] 
          out[[sens]][[sim]][["ICUR"]][arm_i]         <-   out[[sens]][[sim]][["dcosts"]][arm_i]             / out[[sens]][[sim]][["dqalys"]][arm_i] 
          out[[sens]][[sim]][["ICER_undisc"]][arm_i]  <-   out[[sens]][[sim]][["dcosts_undisc"]][arm_i]      / out[[sens]][[sim]][["dlys_undisc"]][arm_i] 
          out[[sens]][[sim]][["ICUR_undisc"]][arm_i]  <-   out[[sens]][[sim]][["dcosts_undisc"]][arm_i]      / out[[sens]][[sim]][["dqalys_undisc"]][arm_i] 
          out[[sens]][[sim]][["INMB"]][arm_i]         <-  out[[sens]][[sim]][["dqalys"]][arm_i] * wtp        - out[[sens]][[sim]][["dcosts"]][arm_i]  
          out[[sens]][[sim]][["INMB_undisc"]][arm_i]  <-  out[[sens]][[sim]][["dqalys_undisc"]][arm_i] * wtp - out[[sens]][[sim]][["dcosts_undisc"]][arm_i]  
          
        } else{
          out[[sens]][[sim]][["ICER"]][arm_i] <- NA
          out[[sens]][[sim]][["ICUR"]][arm_i] <- NA
          out[[sens]][[sim]][["ICER_undisc"]][arm_i]  <- NA
          out[[sens]][[sim]][["ICUR_undisc"]][arm_i]  <- NA
          out[[sens]][[sim]][["INMB"]][arm_i]  <- NA
          out[[sens]][[sim]][["INMB_undisc"]][arm_i]  <- NA
        }
      }
    }
 
  
    data <- data.frame(
      arm           = out[[1]][[1]]$arm_list,
      analysis      = sens,
      analysis_name = out[[sens]][[1]][["sensitivity_name"]],
      costs         = interval_out(out[[sens]],"total_costs",0),
      dcosts        = interval_out(out[[sens]],"dcosts",0),
      lys           = interval_out(out[[sens]],"total_lys",2),
      dlys          = interval_out(out[[sens]],"dlys",3),
      qalys         = interval_out(out[[sens]],"total_qalys",2),
      dqalys        = interval_out(out[[sens]],"dqalys",3),
      ICER          = interval_out(out[[sens]],"ICER",0),
      ICUR          = interval_out(out[[sens]],"ICUR",0),
      INMB          = interval_out(out[[sens]],"INMB",0),
      costs_undisc  = interval_out(out[[sens]],"total_costs_undisc",0),
      dcosts_undisc = interval_out(out[[sens]],"dcosts_undisc",0),
      lys_undisc    = interval_out(out[[sens]],"total_lys_undisc",2),
      dlys_undisc   = interval_out(out[[sens]],"dlys_undisc",3),
      qalys_undisc  = interval_out(out[[sens]],"total_qalys_undisc",2),
      dqalys_undisc = interval_out(out[[sens]],"dqalys_undisc",3),
      ICER_undisc   = interval_out(out[[sens]],"ICER_undisc",0),
      ICUR_undisc   = interval_out(out[[sens]],"ICUR_undisc",0),
      INMB_undisc   = interval_out(out[[sens]],"INMB_undisc",0)
    )
    for (output_i in other_outputs) {
      data[[output_i]] <- interval_out(out[[sens]],output_i,2)
      data[[paste0("d",output_i)]] <-  interval_out(out[[sens]],paste0("d",output_i),3)
    }
    
    names <-  data[,1]
    
    # # Transpose everything other than the first column
    # data <- as.data.frame(as.matrix(t(data[,-1])))
    # 
    # # Assign first column as the column names of the transposed dataframe
    # colnames(data) <- names
    
    data_final <- rbind(data_final,data)
  
  }
  
  data_final <- melt(as.data.table(data_final), id.vars = c("arm", "analysis","analysis_name"))
  
  return(data_final)
  
}


# Extract all specific PSA result -------------------------------------------------------------------------------------------------------------------------

#' Extract PSA results from a treatment
#'
#' @param x The output_sim data frame from the list object returned by `run_sim()`
#' @param element Variable for which PSA results are being extracted (single string)
#'
#' @return A dataframe with PSA results from the specified intervention
#' @export
#'
#' @examples
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' 
#' extract_psa_result(res[[1]],"total_costs")

extract_psa_result <- function(x, element) {
  out <- as.data.frame(do.call(rbind, map(x,element)))
  out$simulation <- 1:nrow(out)
  out$element <- element
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
#' 
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' 
#' @export
#'
#' @examples
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' ceac_des(seq(from=10000,to=500000,by=10000),res)
#' 

ceac_des <- function(wtp, results, interventions = NULL, sensitivity_used = 1) {

  if (is.null(interventions)) {
    interventions <- results[[sensitivity_used]][[1]]$arm_list
  }

  nmb <- data.frame()
  for (comparator in interventions) {

     nmb_i <- t(as.matrix(sapply(wtp, function(wtp_i) wtp_i * extract_psa_result(results[[sensitivity_used]],"total_qalys")[,comparator] - extract_psa_result(results[[sensitivity_used]],"total_costs")[,comparator])))

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
    dplyr::group_by(wtp,iteration) %>%
    dplyr::mutate(best_nmb= ifelse(max(nmb)>0,comparator[nmb==max(nmb)],NA))

  ceac <- nmb %>%
    dplyr::group_by(wtp,comparator) %>%
    dplyr::summarise(prob_best= sum(best_nmb==comparator)/dplyr::n()) %>%
    dplyr::mutate(prob_best = ifelse(is.na(prob_best),0,prob_best))


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
#' 
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' 
#' @export
#'
#' @examples
#' res <- list(list(list(sensitivity_name = "", arm_list = c("int", "noint"
#' ), total_lys = c(int = 9.04687362556945, noint = 9.04687362556945
#' ), total_qalys = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), total_costs = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), total_lys_undisc = c(int = 10.8986618377039, noint = 10.8986618377039
#' ), total_qalys_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), total_costs_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), c_default = c(int = 49921.6357486899, noint = 41225.2544659378
#' ), c_default_undisc = c(int = 59831.3573929783, noint = 49293.1025437205
#' ), q_default = c(int = 6.20743830697466, noint = 6.18115138126336
#' ), q_default_undisc = c(int = 7.50117621700097, noint = 7.47414569286751
#' ), merged_df = list(simulation = 1L, sensitivity = 1L))))
#' 
#' evpi_des(seq(from=10000,to=500000,by=10000),res)
#' 

evpi_des <- function(wtp, results, interventions = NULL, sensitivity_used = 1) {

  if (is.null(interventions)) {
    interventions <- results[[sensitivity_used]][[1]]$arm_list
  }


  nmb <- data.frame()
  for (comparator in interventions) {

    nmb_i <- t(as.matrix(sapply(wtp, function(wtp_i) wtp_i * extract_psa_result(results[[sensitivity_used]],"total_qalys")[,comparator] - extract_psa_result(results[[sensitivity_used]],"total_costs")[,comparator])))

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
    dplyr::group_by(wtp,comparator) %>%
    dplyr::mutate(mean_nmb=mean(nmb)) %>%
    dplyr::group_by(wtp,iteration) %>%
    dplyr::mutate(max_nmb=max(nmb)) %>%
    dplyr::group_by(wtp) %>%
    dplyr::mutate(max_mean_nmb = max(mean_nmb),
           mean_max_nmb = mean(max_nmb)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(evpi = mean_max_nmb - max_mean_nmb) %>%
    dplyr::select(wtp,evpi)


  return(nmb)
}
