# Extracts items and events by looking into assignments, modify_event and new_event

Extracts items and events by looking into assignments, modify_event and
new_event

## Usage

``` r
extract_elements_from_list(node, conditional_flag = FALSE)
```

## Arguments

- node:

  Relevant node within the nested AST list

- conditional_flag:

  Boolean whether the statement is contained within a conditional
  statement

## Value

A data.frame with the relevant item/event, the event where it's
assigned, and whether it's contained within a conditional statement

## Examples

``` r
expr <- substitute({
  
  a <- sum(5+7)
  
  ggplot()
  
  data.frame(x=1,b=2)
  
  list(b=5)
  
  a <- list(s=7)
  
  
  j <- 6
  if(TRUE){modify_event(list(j=5))}
  
  l <- 9
  
  afsa=ifelse(TRUE,"asda",NULL)
  
  
  o_exn = o_exn + 1
  
  a = NULL
  
  b = if(a){"CZ"}else{"AW"}
  
  rnd_prob_exn_sev = runif(1)
  
  exn_sev = rnd_prob_exn_sev <= p_sev
  
  o_exn_mod = o_exn_mod + if(exn_sev) { 0 } else { 1 }
  
  o_exn_sev = o_exn_sev + if(exn_sev) { 1 } else { 0 }
  
  o_rec_time_without_exn = (o_exn == 0) * 1
  
  o_rec_time_without_exn_sev = (o_exn_sev == 0) * 1
  
  o_c_exn = if(exn_sev) { c_sev } else { c_mod }
  
  o_other_c_exn_mod = if(exn_sev) { 0 } else { c_mod }
  
  o_other_c_exn_sev = if(exn_sev) { c_sev } else { 0 }
  
  o_qloss_exn = -if(exn_sev) { q_sev } else { q_mod }
  
  o_other_qloss_exn_mod = -if(exn_sev) { 0 } else { q_mod }
  
  o_other_qloss_exn_sev = -if(exn_sev) { q_sev } else { 0 }
  
  o_qloss_cg_exn = -if(exn_sev) { q_cg_sev } else { q_cg_mod }
  
  o_other_qloss_cg_exn_mod = -if(exn_sev) { 0 } else { q_cg_mod }
  
  o_other_qloss_cg_exn_sev = -if(exn_sev) { q_cg_sev } else { 0 }
  
  o_q = utility
  
  o_other_q_gold1 = if(gold == 1) { utility } else { 0 }
  
  o_other_q_gold2 = if(gold == 2) { utility } else { 0 }
  
  o_other_q_gold3 = if(gold == 3) { utility } else { 0 }
  
  o_other_q_gold4 = if(gold == 4) { utility } else { 0 }
  
  o_other_q_on_dup = if(on_dup) { utility } else { 0 }
  
  n_exn = n_exn + 1
  
  n_exn_mod = n_exn_mod + (1 - exn_sev)
  
  n_exn_sev = n_exn_sev + exn_sev
  
  u_adj_exn_lt = u_adj_exn_lt + if(exn_sev) { u_adj_sev_lt } else { u_adj_mod_lt }
  
  utility = u_gold - u_adj_exn_lt - u_mace_lt
  
  o_rec_utility = utility
  
  rnd_exn = runif(1)
  
  
  if(a==1){
    a=list(6+b)
    
    modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
  } else{
    modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
    if(a>6){
      a=8
    }
    
  }
  
  
  if (sel_resp_incl == 1 & on_dup == 1) {
    
    modify_event(list(e_response = curtime, z = 6))
    
  }
  
})


out <- ast_as_list(expr)

results <- extract_elements_from_list(out)

```
