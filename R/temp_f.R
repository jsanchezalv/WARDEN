library(lobstr)

library(rlang)

getAST <- function(ee) purrr::map_if(as.list(ee), is.call, getAST)

extract_elements_from_list <- function(node, conditional_flag = FALSE) {
  results <- data.frame(
    name = character(),
    definition = character(),
    type = character(),
    conditional_flag = logical(),
    stringsAsFactors = FALSE
  )
  
  # Base case: check if the node is a relevant function
  if (is.list(node)) {
    # Safely extract and convert node[[1]] to character
    func_name <- if (!is.null(node[[1]])) as.character(node[[1]]) else NULL
    
    if (!is.null(func_name) && func_name %in% c("modify_item_seq", "modify_item", "modify_event", "new_event")) {
      # Determine type
      type <- if (func_name %in% c("modify_item_seq", "modify_item")) "item" else "event"
      
      # Extract list elements
      list_expr <- node[[2]]
      if (is.list(list_expr)) {
        for (name in names(list_expr)) {
          definition <- deparse(list_expr[[name]], width.cutoff = 500)
          results <- rbind(results, data.frame(
            name = name,
            definition = definition,
            type = type,
            conditional_flag = conditional_flag,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Check if the node is an `if` block
    if (!is.null(func_name) && func_name == "if") {
      # Process the condition
      conditional_flag <- TRUE
    }
  }
  
  # Recursive case: iterate through list elements
  if (is.list(node)) {
    for (child in node) {
      results <- rbind(results, extract_elements_from_list(child, conditional_flag))
    }
  }
  
  results <- results[!(is.na(results$name) | results$name == "" | is.na(results$definition)), ]
  
  
  return(results)
}



# Function to clean and improve the results
clean_results <- function(results) {
  # Filter out rows with NULL or blank names/definitions
  results <- results[!(is.na(results$name) | results$name == ""), ]
  
  # Replace 'list(...)' in definitions with readable expressions
  results$definition <- sapply(results$definition, function(def) {
    tryCatch({
      # Parse the expression into an R object
      expr <- parse(text = def)[[1]]
      # Use deparse_list for better formatting
      a <- deparse_list(expr)
    }, error = function(e) def) # Fall back to original definition on error
  })
  
  return(results)
}






expr_from_list <- function(lst) {
  if(is.null(lst)){
    return(lst)
  } else if (is.atomic(lst)) {
    return(lst)
  } else if (is.symbol(lst)) {
    return(as.name(lst))
  } else if (is.list(lst) && length(lst) == 1) {
    return(expr_from_list(lst[[1]]))
  } else {
    func_name <- as.character(expr_from_list(lst[[1]]))
    
    # Ensure func_name is a character string
    if (!is.character(func_name)) {
      stop("Function name must be a character string: ", func_name)
    }
    
    args <- lapply(lst[-1], expr_from_list)
    
    # Handle special cases (e.g., operators)
      return(do.call(call, c(func_name, args),quote = TRUE))
  }
}

clean_output <- function(called){
  gsub("    ","",
       gsub("\"","'",
            paste0(
              deparse(called
              ),collapse="")
       )
  )
}


extract_defs <- function(y){
  if(purrr::pluck_depth(y)<=2 & is.null(names(y))){
    NULL
  }else{
    if(length(y[!is.null(names(y))&names(y)!=""])==0){
      lapply(y, extract_defs)
    } else{
      lapply(y[!is.null(names(y)) &names(y)!=""], function(x){
          if(is.null(x)){"NULL"}else{
            clean_output(expr_from_list(x))              
            
          }
      }
    )
    }
  }
}



expr <- substitute({
  
  j <- 6
  if(TRUE){modify_event(list(j=5))}
  
  l <- 9
  
  modify_item(list(afsa=ifelse(TRUE,"asda",NULL)))
  
  modify_item_seq(list(
    
    o_exn = o_exn + 1,
    
    a = NULL,
    
    b = if(a){"CZ"}else{"AW"},
    
    rnd_prob_exn_sev = runif_stream(1, substream_prob_exn_sev),
    
    exn_sev = rnd_prob_exn_sev <= p_sev,
    
    o_exn_mod = o_exn_mod + if(exn_sev) { 0 } else { 1 },
    
    o_exn_sev = o_exn_sev + if(exn_sev) { 1 } else { 0 },
    
    o_rec_time_without_exn = (o_exn == 0) * 1,
    
    o_rec_time_without_exn_sev = (o_exn_sev == 0) * 1,
    
    o_c_exn = if(exn_sev) { c_sev } else { c_mod },
    
    o_other_c_exn_mod = if(exn_sev) { 0 } else { c_mod },
    
    o_other_c_exn_sev = if(exn_sev) { c_sev } else { 0 },
    
    o_qloss_exn = -if(exn_sev) { q_sev } else { q_mod },
    
    o_other_qloss_exn_mod = -if(exn_sev) { 0 } else { q_mod },
    
    o_other_qloss_exn_sev = -if(exn_sev) { q_sev } else { 0 },
    
    o_qloss_cg_exn = -if(exn_sev) { q_cg_sev } else { q_cg_mod },
    
    o_other_qloss_cg_exn_mod = -if(exn_sev) { 0 } else { q_cg_mod },
    
    o_other_qloss_cg_exn_sev = -if(exn_sev) { q_cg_sev } else { 0 },
    
    o_q = utility,
    
    o_other_q_gold1 = if(gold == 1) { utility } else { 0 },
    
    o_other_q_gold2 = if(gold == 2) { utility } else { 0 },
    
    o_other_q_gold3 = if(gold == 3) { utility } else { 0 },
    
    o_other_q_gold4 = if(gold == 4) { utility } else { 0 },
    
    o_other_q_on_dup = if(on_dup) { utility } else { 0 },
    
    n_exn = n_exn + 1,
    
    n_exn_mod = n_exn_mod + (1 - exn_sev),
    
    n_exn_sev = n_exn_sev + exn_sev,
    
    u_adj_exn_lt = u_adj_exn_lt + if(exn_sev) { u_adj_sev_lt } else { u_adj_mod_lt },
    
    utility = u_gold - u_adj_exn_lt - u_mace_lt,
    
    o_rec_utility = utility,
    
    rnd_exn = runif_stream(1, substream_exn)
    
  ))
  
  if(a==1){
    modify_item(list(a=list(6+b)))
    
    modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
  } else{
    modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
    if(a>6){
      modify_item(list(a=8))
    }
    
  }
  
  
  if (sel_resp_incl == 1 & on_dup == 1) {
    
    modify_event(list(e_response = curtime, z = 6))
    
  }
  
})


out <- getAST(expr)

results <- extract_elements_from_list(out)

cleaned_results <- clean_results(results)


a <- unlist(extract_defs(out),recursive=TRUE)
a <- a[!names(a)==""]

results$definition <- a