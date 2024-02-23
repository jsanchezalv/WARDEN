#' Run simulations in parallel mode (at the simulation level)
#'
#' @param arm_list A vector of the names of the interventions evaluated in the simulation
#' @param sensitivity_inputs A list of sensitivity inputs that do not change within a sensitivity in a similar fashion to common_all_inputs, etc
#' @param common_all_inputs A list of inputs common across patients that do not change within a simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param init_event_list A list of initial events and event times. If no initial events are given, a "Start" event at time 0 is created automatically
#' @param evt_react_list A list of event reactions
#' @param util_ongoing_list A list of utilities that are accrued at an ongoing basis
#' @param util_instant_list A list of utilities that are accrued instantaneously at an event
#' @param util_cycle_list A list of utilities that are accrued in cycles
#' @param cost_ongoing_list A list of costs that are accrued at an ongoing basis
#' @param cost_instant_list A list of costs that are accrued instantaneously at an event
#' @param cost_cycle_list A list of costs that are accrued in cycles
#' @param npats The number of patients to be simulated (it will simulate npats * length(arm_list))
#' @param n_sim The number of simulations to run per sensitivity
#' @param psa_bool A boolean to determine if PSA should be conducted. If n_sim > 1 and psa_bool = FALSE, the differences between simulations will be due to sampling
#' @param sensitivity_bool A boolean to determine if Scenarios/DSA should be conducted. 
#' @param sensitivity_names A vector of scenario/DSA names that can be used to select the right sensitivity (e.g., c("Scenario_1", "Scenario_2")). The parameter "sens_name_used" is created from it which corresponds to the one being used for each iteration.
#' @param n_sensitivity Number of sensitivity analysis (DSA or Scenarios) to run. It will be interacted with sensitivity_names argument if not null (n_sensitivityitivity = n_sensitivity * length(sensitivity_names)). For DSA, it should be as many parameters as there are. For scenario, it should be 1.
#' @param ncores The number of cores to use for parallel computing
#' @param drc The discount rate for costs
#' @param drq The discount rate for LYs/QALYs
#' @param input_out A vector of variables to be returned in the output data frame
#' @param ipd A boolean to determine if individual patient data should be returned. If set to false, only the main aggregated outputs will be returned (slightly speeds up code)
#'
#' @return A list of lists with the analysis results
#' @importFrom doParallel registerDoParallel
#'
#' @export
#' 
#' @details This function is slightly different from `run_sim`.
#' `run_sim` allows to run single-core.
#' `run_sim_parallel` allows to use multiple-core at the simulation level,
#' making it more efficient for a large number of simulations relative to `run_sim` (e.g., for  PSA).
#' A list of protected objects that should not be used by the user as input names to avoid the risk of overwriting them is as follows:
#' c("arm", "arm_list", "categories_for_export", "cur_evtlist", "curtime", "evt", "i", "prevtime", "sens", "simulation", "sens_name_used","list_env","uc_lists") 
#' 
#' @examples
#' \dontrun{
#' run_sim_parallel(arm_list=c("int","noint"),
#' common_all_inputs = common_all_inputs,
#' common_pt_inputs = common_pt_inputs,
#' unique_pt_inputs = unique_pt_inputs,
#' init_event_list = init_event_list,
#' evt_react_list = evt_react_list,
#' util_ongoing_list = util_ongoing_list,
#' util_instant_list = util_instant_list,
#' cost_ongoing_list = cost_ongoing_list,
#' cost_instant_list = cost_instant_list,
#' npats = 500,
#' n_sim = 1,
#' psa_bool = FALSE,
#' ncores = 3,
#' drc = 0.035,
#' drq = 0.035,
#' ipd = TRUE)
#' }

run_sim_parallel <- function(arm_list=c("int","noint"),
                   sensitivity_inputs=NULL,
                   common_all_inputs=NULL,
                   common_pt_inputs=NULL,
                   unique_pt_inputs=NULL,
                   init_event_list = NULL,
                   evt_react_list = evt_react_list,
                   util_ongoing_list = NULL,
                   util_instant_list = NULL,
                   util_cycle_list = NULL,
                   cost_ongoing_list = NULL,
                   cost_instant_list = NULL,
                   cost_cycle_list = NULL,
                   npats=500,
                   n_sim=1,
                   psa_bool = NULL,
                   sensitivity_bool = FALSE,
                   sensitivity_names = NULL,
                   n_sensitivity = 1,
                   ncores=1,
                   drc=0.035,
                   drq=0.035,
                   input_out = NULL,
                   ipd = TRUE){


# Set-up basics -----------------------------------------------------------

  cl <- parallel::makeCluster(ncores)
  registerDoParallel(cl)
  
  arm_list <- arm_list #this is done as otherwise there are problems passing arguments from function to function
  
  #get cost/utility categories
  categories_costs_ongoing <- unlist(unique(lapply(names(cost_ongoing_list), function(n) cost_ongoing_list[[n]][["category"]])))
  categories_costs_instant <- unlist(unique(lapply(names(cost_instant_list), function(n) cost_instant_list[[n]][["category"]])))
  categories_costs_cycle   <- unlist(unique(lapply(names(cost_cycle_list), function(n) cost_cycle_list[[n]][["category"]])))

  categories_utilities_ongoing <- unlist(unique(lapply(names(util_ongoing_list), function(n) util_ongoing_list[[n]][["category"]])))
  categories_utilities_instant <- unlist(unique(lapply(names(util_instant_list), function(n) util_instant_list[[n]][["category"]])))
  categories_utilities_cycle   <- unlist(unique(lapply(names(util_cycle_list), function(n) util_cycle_list[[n]][["category"]])))
  
  categories_for_export <-c("categories_costs_ongoing",
                            "categories_costs_instant",
                            "categories_costs_cycle",
                            "categories_utilities_ongoing",
                            "categories_utilities_instant",
                            "categories_utilities_cycle"
                            )
  
  #Remove NULL values
  categories_for_export <- c(if(!is.null(categories_costs_ongoing)){paste(categories_costs_ongoing,c("ongoing","ongoing_undisc"),sep="_")},
                             if(!is.null(categories_costs_instant)){paste(categories_costs_instant,c("instant","instant_undisc"),sep="_")},
                             if(!is.null(categories_costs_cycle)){paste(categories_costs_cycle,c("cycle","cycle_undisc"),sep="_")},
                             if(!is.null(categories_utilities_ongoing)){paste(categories_utilities_ongoing,c("ongoing","ongoing_undisc"),sep="_")},
                             if(!is.null(categories_utilities_instant)){paste(categories_utilities_instant,c("instant","instant_undisc"),sep="_")},
                             if(!is.null(categories_utilities_cycle)){paste(categories_utilities_cycle,c("cycle","cycle_undisc"),sep="_")}
                            )
  output_sim <- list()

  start_time_simulations <-  proc.time()
  
# Sensitivity loop ---------------------------------------------------------
  
  if (is.null(sensitivity_names)) {
    length_sensitivities <- n_sensitivity
  } else{
    length_sensitivities <- n_sensitivity * length(sensitivity_names)
  }
  
  for (sens in 1:length_sensitivities) {
    print(paste0("Sensitivity number: ",sens))
    start_time <-  proc.time()

    output_sim[[sens]] <- list() #initialize sensitivity lists
    
    #e.g., if length_sensitivities is 50 (25 param x 2 DSAs) then take at each iteration the divisor to see which name should be applied
    if (!is.null(sensitivity_names)) {
    sens_name_used <- sensitivity_names[ceiling(sens/n_sensitivity)] 
    } else{
    sens_name_used <- NULL
    }
    
    
    input_list_sens <- list(drc = drc,
                       drq = drq,
                       psa_bool = psa_bool,
                       init_event_list = init_event_list,
                       evt_react_list = evt_react_list,
                       uc_lists = list(util_ongoing_list = util_ongoing_list,
                                       util_instant_list = util_instant_list,
                                       util_cycle_list = util_cycle_list,
                                       cost_ongoing_list = cost_ongoing_list,
                                       cost_instant_list = cost_instant_list,
                                       cost_cycle_list = cost_cycle_list,
                                       cost_categories_ongoing = categories_costs_ongoing,
                                       l_cost_categories_ongoing = length(categories_costs_ongoing),
                                       cost_categories_instant = categories_costs_instant,
                                       l_cost_categories_instant = length(categories_costs_instant),
                                       cost_categories_cycle = categories_costs_cycle,
                                       l_cost_categories_cycle = length(categories_costs_cycle),
                                       util_categories_ongoing = categories_utilities_ongoing,
                                       l_util_categories_ongoing = length(categories_utilities_ongoing),
                                       util_categories_instant = categories_utilities_instant,
                                       l_util_categories_instant = length(categories_utilities_instant),
                                       util_categories_cycle = categories_utilities_cycle,
                                       l_util_categories_cycle = length(categories_utilities_cycle)
                       ),
                       input_out = unique(c(input_out,categories_for_export)),
                       categories_for_export = categories_for_export,
                       ipd = ipd,
                       arm_list = arm_list,
                       npats = npats,
                       n_sim = n_sim,
                       n_sensitivity = n_sensitivity,
                       sens = sens,
                       sensitivity_names = sensitivity_names,
                       sens_name_used = sens_name_used
                      )
    
    # Draw Common parameters  -------------------------------
    if(!is.null(sensitivity_inputs)){
      for (inp in 1:length(sensitivity_inputs)) {
        set.seed(sens)
        list.sensitivity_inputs <- lapply(sensitivity_inputs[inp],function(x) eval(x, input_list_sens))
        if ((!is.null(names(list.sensitivity_inputs[[1]]))) & sens==1) {
          warning("Item ", names(list.sensitivity_inputs), " is named. It is strongly advised to assign unnamed objects if they are going to be processed in the model, as they could generate errors.")
        }
        input_list_sens <- c(input_list_sens,list.sensitivity_inputs)
      }
    }

# Simulation loop ---------------------------------------------------------

    output_sim[[sens]] <- vector("list", length=n_sim) # empty list with n_sim elements
    # Outer loop, repeat for each patient
    output_sim[[sens]] <- foreach(simulation = 1:n_sim,
                         .packages = (.packages()),
                         .export = unique(c("input_list_sens",ls(.GlobalEnv),ls(parent.env(environment())),ls(environment()))),
                         .combine = 'c') %dopar% {

      print(paste0("Simulation number: ",simulation))
      input_list <- c(input_list_sens,
                      simulation = simulation)
  
      # Draw Common parameters  -------------------------------
      if(!is.null(common_all_inputs)){
        for (inp in 1:length(common_all_inputs)) {
          set.seed(simulation)
          list.common_all_inputs <- lapply(common_all_inputs[inp],function(x) eval(x, input_list))
          if ((!is.null(names(list.common_all_inputs[[1]]))) & simulation==1 & sens==1) {
            warning("Item ", names(list.common_all_inputs), " is named. It is strongly advised to assign unnamed objects if they are going to be processed in the model, as they could generate errors.")
          }
          input_list <- c(input_list,list.common_all_inputs)
        }
      }
  
      #Make sure there are no duplicated inputs in the model, if so, take the last one
      duplic <- duplicated(names(input_list),fromLast = T)
      if (sum(duplic)>0 & simulation==1 & sens==1) { warning("Duplicated items detected, using the last one added")  }
      input_list <- input_list[!duplic]
  
      # Run engine ----------------------------------------------------------
  
        final_output <- run_engine(arm_list=arm_list,
                                        common_pt_inputs=common_pt_inputs,
                                        unique_pt_inputs=unique_pt_inputs,
                                        input_list = input_list)                    # run simulation
      
  
      if (input_list$ipd==TRUE) {
  
        final_output$merged_df$simulation <- simulation
        final_output$merged_df$sensitivity <- sens
      }
  
      return(list(final_output))
      
      print(paste0("Time to run simulation ", simulation,": ",  round(proc.time()[3]- start_time[3] , 2 ), "s"))
    }
    
    print(paste0("Time to run analysis ", sens,": ",  round(proc.time()[3]- start_time[3] , 2 ), "s"))
    
  }
  print(paste0("Total time to run: ",  round(proc.time()[3]- start_time_simulations[3] , 2), "s"))
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals) 
  parallel::stopCluster(cl)
  

  # Export results ----------------------------------------------------------


  results <- output_sim

  return(results)

}
