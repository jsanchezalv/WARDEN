
# Global variables for CRAN check -----------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      c('simulation')
    )) 
}

#' Run simulations in parallel mode (at the simulation level)
#'
#' @param arm_list A vector of the names of the interventions evaluated in the simulation
#' @param sensitivity_inputs A list of sensitivity inputs that do not change within a sensitivity in a similar fashion to common_all_inputs, etc
#' @param common_all_inputs A list of inputs common across patients that do not change within a simulation
#' @param common_pt_inputs A list of inputs that change across patients but are not affected by the intervention
#' @param unique_pt_inputs A list of inputs that change across each intervention
#' @param init_event_list A list of initial events and event times. If no initial events are given, a "Start" event at time 0 is created automatically
#' @param evt_react_list A list of event reactions
#' @param util_ongoing_list Vector of QALY named variables that are accrued at an ongoing basis (discounted using drq)
#' @param util_instant_list Vector of QALY named variables that are accrued instantaneously at an event (discounted using drq)
#' @param util_cycle_list Vector of QALY named variables that are accrued in cycles (discounted using drq)
#' @param cost_ongoing_list Vector of cost named variables that are accrued at an ongoing basis (discounted using drc)
#' @param cost_instant_list Vector of cost named variables that are accrued instantaneously at an event (discounted using drc)
#' @param cost_cycle_list Vector of cost named variables that are accrued in cycles (discounted using drc)
#' @param other_ongoing_list Vector of other named variables that are accrued at an ongoing basis (discounted using drq)
#' @param other_instant_list Vector of other named variables that are accrued instantaneously at an event (discounted using drq)
#' @param npats The number of patients to be simulated (it will simulate npats * length(arm_list))
#' @param n_sim The number of simulations to run per sensitivity
#' @param psa_bool A boolean to determine if PSA should be conducted. If n_sim > 1 and psa_bool = FALSE, the differences between simulations will be due to sampling
#' @param sensitivity_bool A boolean to determine if Scenarios/DSA should be conducted. 
#' @param sensitivity_names A vector of scenario/DSA names that can be used to select the right sensitivity (e.g., c("Scenario_1", "Scenario_2")). The parameter "sens_name_used" is created from it which corresponds to the one being used for each iteration.
#' @param n_sensitivity Number of sensitivity analysis (DSA or Scenarios) to run. It will be interacted with sensitivity_names argument if not null (n_sensitivityitivity = n_sensitivity * length(sensitivity_names)). For DSA, it should be as many parameters as there are. For scenario, it should be 1.
#' @param ncores The number of cores to use for parallel computing
#' @param input_out A vector of variables to be returned in the output data frame
#' @param ipd Integer taking value 0 if no IPD data returned, 1 for full IPD data returned, and 2 IPD data but aggregating events
#' @param constrained Boolean, FALSE by default, which runs the simulation with patients not interacting with each other, TRUE if resources are shared within an arm (allows constrained resources)
#' @param timed_freq If NULL, it does not produce any timed outputs. Otherwise should be a number (e.g., every 1 year)
#' @param debug If TRUE, will generate a log file
#' @param accum_backwards If TRUE, the ongoing accumulators will count backwards (i.e., the current value is applied until the previous update). If FALSE, the current value is applied between the current event and the next time it is updated. If TRUE, user must use `modify_item` and `modify_item_seq` or results will be incorrect.
#' @param continue_on_error If TRUE, on error  at patient stage will attempt to continue to the next simulation (only works if n_sim and/or n_sensitivity are > 1, not at the patient level)
#' @param seed Starting seed to be used for the whole analysis. If null, it's set to 1 by default.
#'
#' @return A list of lists with the analysis results
#' @importFrom doFuture `%dofuture%`
#' @importFrom future plan
#' @importFrom future multisession
#' @importFrom foreach foreach
#' @importFrom progressr with_progress
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom progressr handler_txtprogressbar
#' @importFrom rlang env_clone
#'
#' @export
#' 
#' @details This function is slightly different from `run_sim`.
#' `run_sim` allows to run single-core.
#' `run_sim_parallel` allows to use multiple-core at the simulation level,
#' making it more efficient for a large number of simulations relative to `run_sim` (e.g., for  PSA).
#' 
#' Event ties are processed in the order declared within the `init_event_list` argument (`evts` argument within the first sublist of that object).
#' To do so, the program automatically adds a sequence from to 0 to the (number of events - 1) times 1e-10 to add to the event times when selecting the event with minimum time.
#' This time has been selected as it's relatively small yet not so small as to be ignored by which.min (see .Machine for more details)
#'  
#' A list of protected objects that should not be used by the user as input names or in the global environment to avoid the risk of overwriting them is as follows:
#' c("arm", "arm_list", "categories_for_export", "cur_evtlist", "curtime", "evt", "i", "prevtime", "sens", "simulation", "sens_name_used","list_env","uc_lists","npats","ipd").
#' 
#' The engine uses the L'Ecuyer-CMRG for the random number generator.
#' Note that if ncores > 1, then results per simulation will only be exactly replicable if using run_sim_parallel 
#' (as seeds are automatically transformed to be seven integer seeds -i.e, L'Ecuyer-CMRG seeds-)
#' Note that the random seeds are set to be unique in their category (i.e., at patient level, patient-arm level, etc.)
#' 
#'  If no `drc` or `drq` parameters are passed within `sensitivity` or `common_all` input lists, these are assigned a default value 0.03 for discounting costs, QALYs and others.
#'
#' Ongoing items will look backward to the last time updated when performing the discounting and accumulation. 
#' This means that the user does not necessarily need to keep updating the value, but only add it when the value 
#' changes looking forward (e.g., o_q = utility at event 1, at event 2 utility does not change, but at event 3 it does, 
#' so we want to make sure to add o_q = utility at event 3 before updating utility. The program will automatically 
#' look back until event 1). Note that in previous versions of the package backward was the default, and now this has switched to forward.
#' 
#' If using `accum_backwards = TRUE`, then it is mandatory for the user to use `modify_item` and `modify_item_seq` in event reactions,
#'   as the standard assignment approach (e.g., `a <- 5`) will not calculate the right results, particularly in the presence of
#'   conditional statements.
#' 
#'  If the `cycle` lists are used, then it is expected the user will declare as well the name of the variable
#'   pasted with `cycle_l` and `cycle_starttime` (e.g., c_default_cycle_l and c_default_cycle_starttime) to 
#'   ensure the discounting can be computed using cycles, with cycle_l being the cycle length, and cycle_starttime 
#'   being the starting time in which the variable started counting. Optionally,  `max_cycles` must also be added (if no
#'   maximum number of cycles, it should be set equal to NA).
#'     
#'  `debug = TRUE` will export a log file with the timestamp up the error in the main working directory. Note that
#'  using this mode without modify_item or modify_item_seq may lead to inaccuracies if assignments are done in non-standard ways,
#'  as the AST may not catch all the relevant assignments (e.g., an assigment like assign(paste("x_",i),5)
#'   in a loop will not be identified, unless using modify_item()).
#'   
#'   If `continue_on_error` is set to FALSE, it will only export analysis level inputs due to the parallel engine
#'    (use single-engine for those inputs) `continue_on_error` will skip the current simulation 
#'    (so it won't continue for the rest of patient-arms) if TRUE. 
#'  Note that this will make the progress bar not correct, as a set of patients that were expected to be run is not.
#'
#' @examples
#' library(magrittr)
#' common_all_inputs <-add_item(
#' util.sick = 0.8,
#' util.sicker = 0.5,
#' cost.sick = 3000,
#' cost.sicker = 7000,
#' cost.int = 1000,
#' coef_noint = log(0.2),
#' HR_int = 0.8,
#' drc = 0.035, #different values than what's assumed by default
#' drq = 0.035,
#' random_seed_sicker_i = sample.int(100000,5,replace = FALSE)
#' )
#' 
#' common_pt_inputs <- add_item(death= max(0.0000001,rnorm(n=1, mean=12, sd=3))) 
#' 
#' unique_pt_inputs <- add_item(fl.sick = 1,
#'                              q_default = util.sick,
#'                              c_default = cost.sick + if(arm=="int"){cost.int}else{0}) 
#'                              
#' init_event_list <- 
#' add_tte(arm=c("noint","int"), evts = c("sick","sicker","death") ,input={
#'   sick <- 0
#'   sicker <- draw_tte(1,dist="exp",
#'    coef1=coef_noint, beta_tx = ifelse(arm=="int",HR_int,1),
#'    seed = random_seed_sicker_i[i])
#'   
#' })   
#' 
#' evt_react_list <-
#' add_reactevt(name_evt = "sick",
#'              input = {}) %>%
#'   add_reactevt(name_evt = "sicker",
#'                input = {
#'                  modify_item(list(q_default = util.sicker,
#'                                   c_default = cost.sicker + if(arm=="int"){cost.int}else{0},
#'                                   fl.sick = 0)) 
#'                }) %>%
#'   add_reactevt(name_evt = "death",
#'                input = {
#'                  modify_item(list(q_default = 0,
#'                                   c_default = 0, 
#'                                   curtime = Inf)) 
#'                }) 
#'                
#' util_ongoing <- "q_default"
#' cost_ongoing <- "c_default"
#'                           
#' 
#' run_sim_parallel(arm_list=c("int","noint"),
#' common_all_inputs = common_all_inputs,
#' common_pt_inputs = common_pt_inputs,
#' unique_pt_inputs = unique_pt_inputs,
#' init_event_list = init_event_list,
#' evt_react_list = evt_react_list,
#' util_ongoing_list = util_ongoing,
#' cost_ongoing_list = cost_ongoing,
#' npats = 2,
#' n_sim = 1,
#' psa_bool = FALSE,
#' ipd = 1,
#' ncores = 1)
#' 

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
                             other_ongoing_list = NULL,
                             other_instant_list = NULL,
                             npats=500,
                             n_sim=1,
                             psa_bool = NULL,
                             sensitivity_bool = FALSE,
                             sensitivity_names = NULL,
                             n_sensitivity = 1,
                             ncores=1,
                             input_out = character(),
                             ipd = 1,
                             constrained = FALSE,
                             timed_freq = NULL,
                             debug = FALSE,
                             accum_backwards = FALSE,
                             continue_on_error = FALSE,
                             seed = NULL){
  
  
  # # ---- Error-context beacon (for clean messages; keeps base traceback) ----
  .warden_ctx <- new.env(parent = emptyenv())
  .warden_ctx$last <- NULL
  .set_last_ctx(stage="Error in setup:initial code setup")
  
  log_sink <- new.env(parent = emptyenv())
  log_sink$entries <- list()
  
  on.exit({
    do.call(RNGkind, as.list(rng_kind_store))
    if (debug) export_log(log_sink$entries,paste0("log_model_",format(Sys.time(), "%Y_%m_%d_%Hh_%mm_%Ss"),".txt"))
    ctx <- .warden_ctx$last
    
    message(
      paste0(ctx$stage,"; ",
             ifelse(!is.na(ctx$sens),paste0("Analysis ", ctx$sens,"; "),""),
             ifelse(!is.na(ctx$simulation),paste0("Simulation ", ctx$simulation,"; "),""),
             ifelse(!is.na(ctx$arm),paste0("Arm ", ctx$arm,"; "),""),
             ifelse(!is.na(ctx$patient_id),paste0("Patient ", ctx$patient_id,"; "),""),
             ifelse(!is.na(ctx$event),paste0("Event ", ctx$event,"; "),""),
             ifelse(!is.na(ctx$time),paste0("Time ", ctx$time,"; "),"")
      ))
    
    if(ctx$stage!="Simulation finalized"){message("An error occurred somewhere. Traceback available through traceback()")}
    
  }, add=TRUE)
  
  
  # Set-up basics -----------------------------------------------------------
  
  #Store original rng kind and use L'Ecuyer-CMRG 
  rng_kind_store <- RNGkind()
  RNGkind("L'Ecuyer-CMRG")
  
  
  #Stop simulation if forbidden objects are used.
  list_forbidden_names <-  c("arm", "arm_list", "categories_for_export", "cur_evtlist", "curtime", "evt", "i", "prevtime", "sens", "simulation", "sens_name_used","list_env","uc_lists","npats","ipd") 
  matched_list_forbidden <- which(list_forbidden_names %in% c(
    ls(parent.frame()),
    names(common_all_inputs),
    names(sensitivity_inputs),
    names(common_pt_inputs),
    names(unique_pt_inputs),
    names(arm_list),
    names(evt_react_list)
  )
  )
  
  if(  length(matched_list_forbidden)>0){
    stop(paste0("Name(s) or object: `", list_forbidden_names[matched_list_forbidden],"` defined belong to the list of forbidden names, which can cause issues in the model.
    This is likely caused by an object being declared in the parent/global environment (e.g., `i`) or in the inputs. 
    You can use functions like find('i') to help you locate them. 
    Please remove or rename those objects. See run_sim or run_sim_parallel for a full list of these names.\n"))
  }
  
  
  if(length(names(evt_react_list)) != length(init_event_list[[1]][["evts"]])){
    dif <- setdiff(names(evt_react_list),init_event_list[[1]][["evts"]])
    stop(paste0("Number of events defined in init_event_list first sublist `evts` is not equal to the number of events defined in evt_react_list.
         Please ensure both have equal length, ",dif," is/are missing from one of the lists."))
  }
  
  plan(multisession, workers = ncores)
  
  arm_list <- arm_list #this is done as otherwise there are problems passing arguments from function to function
  
  #get cost/utility categories
  categories_costs_ongoing <- cost_ongoing_list
  categories_costs_instant <- cost_instant_list
  categories_costs_cycle   <- cost_cycle_list
  
  categories_utilities_ongoing <- util_ongoing_list
  categories_utilities_instant <- util_instant_list
  categories_utilities_cycle   <- util_cycle_list
  
  categories_other_ongoing <- other_ongoing_list
  categories_other_instant <- other_instant_list
  
  
  #Remove NULL values
  categories_for_export <- unique(
    c(if(!is.null(categories_costs_ongoing)){categories_costs_ongoing},
      if(!is.null(categories_costs_instant)){categories_costs_instant},
      if(!is.null(categories_costs_cycle)){c(categories_costs_cycle,paste0(categories_costs_cycle,"_cycle_l"),paste0(categories_costs_cycle,"_cycle_starttime"),paste0(categories_costs_cycle,"_max_cycles"))},
      if(!is.null(categories_utilities_ongoing)){categories_utilities_ongoing},
      if(!is.null(categories_utilities_instant)){categories_utilities_instant},
      if(!is.null(categories_utilities_cycle)){c(categories_utilities_cycle,paste0(categories_utilities_cycle,"_cycle_l"),paste0(categories_utilities_cycle,"_cycle_starttime"),paste0(categories_utilities_cycle,"_max_cycles"))},
      if(!is.null(categories_other_ongoing)){categories_other_ongoing},
      if(!is.null(categories_other_instant)){categories_other_instant}
    ))
  
  
  env_setup_sens <- is.language(sensitivity_inputs)
  env_setup_sim <- is.language(common_all_inputs)
  env_setup_pt <- is.language(common_pt_inputs)
  env_setup_arm <- is.language(unique_pt_inputs)

  output_sim <- list()
  
  final_log <- list()
  log_list <- list()
  .skip_to_next <- FALSE
  flag_noerr_sim <- TRUE
  
  start_time <-  proc.time()
  
  # Analysis loop ---------------------------------------------------------
  if (is.null(sensitivity_names)) {
    length_sensitivities <- n_sensitivity
  } else{
    length_sensitivities <- n_sensitivity * length(sensitivity_names)
  }
  
  progressr::handlers(progressr::handler_txtprogressbar(width=100))
  options(progressr.interrupts = FALSE)
  progressr::with_progress({
    pb <- progressr::progressor(min(npats*length_sensitivities*n_sim,50)) 
    
  for (sens in 1:length_sensitivities) {
    .skip_to_next <- FALSE
    
    message(paste0("Analysis number: ",sens))

    start_time_analysis <-  proc.time()
    
    output_sim[[sens]] <- list() #initialize analysis lists
    
    #e.g., if length_sensitivities is 50 (25 param x 2 DSAs) then take at each iteration the divisor to see which name should be applied
    if (!is.null(sensitivity_names)) {
    sens_name_used <- sensitivity_names[ceiling(sens/n_sensitivity)] 
    } else{
    sens_name_used <- ""
    }
    
    if(debug){
      for (evt in names(evt_react_list)) {
        evt_react_list[[evt]]$debug_vars <- extract_assignment_targets(evt_react_list[[evt]]$react)
      }
    }
    
    input_list_sens <- list(
                       psa_bool = psa_bool,
                       init_event_list = init_event_list,
                       precision_times = if(!is.null(init_event_list)){setNames((seq_along(init_event_list[[1]][["evts"]])-1)*1e-10,
                                                                                init_event_list[[1]][["evts"]])},
                       evt_react_list = evt_react_list,
                       uc_lists = list(util_ongoing_list = util_ongoing_list,
                                       util_instant_list = util_instant_list,
                                       util_cycle_list = util_cycle_list,
                                       cost_ongoing_list = cost_ongoing_list,
                                       cost_instant_list = cost_instant_list,
                                       cost_cycle_list = cost_cycle_list,
                                       other_ongoing_list = other_ongoing_list,
                                       other_instant_list = other_instant_list,
                                       ongoing_inputs = c(
                                         if(!is.null(categories_costs_ongoing)){categories_costs_ongoing},
                                         if(!is.null(categories_utilities_ongoing)){categories_utilities_ongoing},
                                         if(!is.null(categories_other_ongoing)){categories_other_ongoing}
                                       ),
                                       instant_inputs = c(
                                         if(!is.null(categories_costs_instant)){categories_costs_instant},
                                         if(!is.null(categories_utilities_instant)){categories_utilities_instant},
                                         if(!is.null(categories_other_instant)){categories_other_instant}
                                       ),
                                       cycle_inputs = c(
                                         if(!is.null(categories_costs_cycle)){categories_costs_cycle},
                                         if(!is.null(categories_utilities_cycle)){categories_utilities_cycle}
                                       ),
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
                                       l_util_categories_cycle = length(categories_utilities_cycle),
                                       other_categories_ongoing = categories_other_ongoing,
                                       l_other_categories_ongoing = length(categories_other_ongoing),
                                       other_categories_instant = categories_other_instant,
                                       l_other_categories_instant = length(categories_other_instant)
                       ),
                       input_out = unique(c(input_out,categories_for_export)),
                       categories_for_export = categories_for_export,
                       ipd = ipd,
                       arm_list = arm_list,
                       npats = npats,
                       n_sim = n_sim,
                       sensitivity_bool = sensitivity_bool,
                       n_sensitivity = n_sensitivity,
                       length_sensitivities = length_sensitivities,
                       sens = sens,
                       sensitivity_names = sensitivity_names,
                       sens_name_used = sens_name_used,
                       timed_freq = timed_freq,
                       debug = debug,
                       accum_backwards = accum_backwards,
                       continue_on_error = continue_on_error,
                       log_list = list(),
                       env_setup_sens = env_setup_sens,
                       env_setup_sim = env_setup_sim,
                       env_setup_pt = env_setup_pt,
                       env_setup_arm = env_setup_arm
                      )
    
    if(is.null(seed)){
      seed <- 1
    }
    set.seed(seed)
    
    # Draw Common parameters  -------------------------------
    input_list_sens <- as.environment(input_list_sens)
    parent.env(input_list_sens) <- environment()
    .set_last_ctx(stage="Error in setup:sensitivity_inputs", sens=sens)
    
    # Draw Common parameters  -------------------------------
    if(!is.null(sensitivity_inputs)){
      
      on_error_check({
        if(env_setup_sens){
          load_inputs2(inputs = input_list_sens,list_uneval_inputs = sensitivity_inputs)
        } else{
          input_list_sens <- as.environment(
            load_inputs(inputs = as.list(input_list_sens),
                        list_uneval_inputs = sensitivity_inputs)
          )
          parent.env(input_list_sens) <- environment()
          
        }
      })
      if(.skip_to_next){next}
      
      if(input_list_sens$debug){ 
        dump_info <- debug_inputs(NULL,input_list_sens)
        
        
        names(dump_info) <- paste0("Analysis: ", input_list_sens$sens," ", input_list_sens$sens_name_used,
                                   "; Structural"
        )
        
        log_list <- c(log_list,dump_info)
        log_add(dump_info)
        
      }
      
    }
    
    
# Simulation loop ---------------------------------------------------------
    output_sim[[sens]] <- vector("list", length=n_sim) # empty list with n_sim elements
    # Outer loop, repeat for each patient
    
    exported_items <- unique(c("input_list_sens",ls(.GlobalEnv),ls(parent.env(environment()), all.names = TRUE),ls(environment(), all.names = TRUE)))
    options(future.rng.onMisuse = "ignore")
    
    output_sim[[sens]] <- foreach(simulation = 1:n_sim,
                         # .options.future = list(seed = TRUE),
                         .options.future = list(packages = .packages()),
                         .options.future = list(globals=structure(TRUE,add = exported_items)),
                         .combine = 'c') %dofuture% {
      .set_last_ctx(stage="Error in setup:simulation_start", sens=sens, simulation=simulation)
                           
      RNGkind("L'Ecuyer-CMRG") #repeat this here so parallel::RngStream does not malfunction
                            
      message(paste0("Simulation number: ",simulation))
      
                           
      start_time_sim <-  proc.time()
                           
                           
      input_list <- rlang::env_clone(input_list_sens , parent.env(input_list_sens))
      input_list$simulation <- simulation
      
      set.seed(simulation*1007*seed)
      .set_last_ctx(stage="setup:common_all_inputs", sens=sens, simulation=simulation)
      
      # Draw Common parameters  -------------------------------
      if(!is.null(common_all_inputs)){
        
        if(env_setup_sim){
          on_error_check({
            if(env_setup_sim){
              load_inputs2(inputs = input_list,list_uneval_inputs = common_all_inputs)
            } else{
              input_list <- as.environment(
                load_inputs(inputs = as.list(input_list),
                            list_uneval_inputs = common_all_inputs)
              )
              parent.env(input_list) <- parent.env(input_list_sens)
            }
          })
          if(.skip_to_next){return(NULL)}
        
        if(input_list_sens$debug){ 
          dump_info <- debug_inputs(input_list_sens,input_list)
          
          
          names(dump_info) <- paste0("Analysis: ", input_list$sens," ", input_list$sens_name_used,
                                     "; Sim: ", input_list$simulation,
                                     "; Statics"
          )
          
          log_list <- c(log_list,dump_info)
          log_add(dump_info)
          
        }
      }
  
      if(is.null(input_list$drc)){input_list$drc <- 0.03}
      if(is.null(input_list$drq)){input_list$drq <- 0.03}
      
      # Run engine ----------------------------------------------------------
        .set_last_ctx(stage="Error in engine:start", sens=sens, simulation=simulation)
        on_error_check({
          if(constrained){
            final_output <- run_engine_constrained(arm_list=arm_list,
                                                   common_pt_inputs=common_pt_inputs,
                                                   unique_pt_inputs=unique_pt_inputs,
                                                   input_list = input_list,
                                                   pb = pb,
                                                   seed = seed,
                                                   .warden_ctx = .warden_ctx)  
          } else{
            final_output <- run_engine(arm_list=arm_list,
                                       common_pt_inputs=common_pt_inputs,
                                       unique_pt_inputs=unique_pt_inputs,
                                       input_list = input_list,
                                       pb = pb,
                                       seed = seed,
                                       .warden_ctx = .warden_ctx)   
          }
        })
        if(.skip_to_next){return(NULL)}
      
      if (input_list$ipd>0) {
  
        final_output$merged_df$simulation <- simulation
        final_output$merged_df$sensitivity <- sens
      }
      
      final_output <- c(list(sensitivity_name = sens_name_used), final_output)
      
      if (debug) final_output$log_list <- log_sink$entries
      
      return(list(final_output))
      
      message(paste0("Time to run simulation ", simulation,": ",  round(proc.time()[3]- start_time_sim[3] , 2 ), "s"))
      

     }
    

    message(paste0("Time to run analysis ", sens,": ",  round(proc.time()[3]- start_time_analysis[3] , 2 ), "s"))
    

  }
  message(paste0("Total time to run: ",  round(proc.time()[3]- start_time[3] , 2), "s"))
 
  
  

  results <- output_sim
  
    }
  }, enable=TRUE, cleanup = TRUE) 
  
  .set_last_ctx(stage="Simulation finalized")
  
  
  return(results)

}
