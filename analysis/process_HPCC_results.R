# process HPCC results

library(tidyverse)
library(dplyr)
library(tidyr)



# load HPCC results
load(here::here("data", "HPCC_results_v2.Rdata"))

HPCC_results <- all_cases |> 
  mutate(log_estimates = as.numeric(log_estimates),
         log_scale_SE = as.numeric(log_scale_SE),
         job_no = as.factor(job_no),
         sim = as.factor(sim))


# compute relative error  --------------------------------------------------------------------


# # add in theta4 and theta5 for each sim for cases that don't estimate it
# thetas_bind <- tibble(parameter = rep(c("theta5", "theta6"), each = 4000),
#                       log_estimates = rep(NA_real_, 8000),
#                       log_scale_SE = rep(NA_real_, each = 8000),
#                       case = rep(c("a", "b", "e", "f"), 2000))
# add_thetas_not_estimated <- pc_19 |> 
#   select(-c(job_no, sim)) |> 
#   bind_rows(thetas_bind) |> arrange(case)
# 
# 
# # this is now the HPCC results table but with theta5 and theta6
# # placeholders for all cases not estimating those correlations - a, b, e, and f
# HPCC_results_new <- bind_rows(pc_23, add_thetas_not_estimated) |> 
#   select(-c(job_no, sim))

# load in true estimates
load(here::here("data", "best_fitting_model_with_PIs.Rdata"))

true_values <- probability_intervals |> 
  mutate(true_values = log_estimates) |> 
  select(parameter, true_values)


relative_error <- HPCC_results |>
  drop_na(parameter) |> 
  left_join(true_values,
            by = "parameter") |>
  mutate(true_values = case_when(case %in% c("a", "c", "e", "g") & parameter == "theta5" ~ 0.001,
                                 case %in% c("a", "c", "e", "g") & parameter == "theta6" ~ 0.001,
                                 case %in% c("b", "d", "f", "h") & parameter == "theta5" ~ 0.58,
                                 case %in% c("b", "d", "f", "h") & parameter == "theta6" ~ 0.58,
                                 case %in% c("a", "c", "e", "g") & parameter == "Linf_L1_pop_cor" ~ 0.001,
                                 case %in% c("a", "c", "e", "g") & parameter == "L1_K_pop_cor" ~ 0.001,
                                 case %in% c("b", "d", "f", "h") & parameter == "Linf_L1_pop_cor" ~ 0.58,
                                 case %in% c("b", "d", "f", "h") & parameter == "L1_K_pop_cor" ~ 0.58,
                                 TRUE ~ true_values),
         relative_error = if_else(is.na(log_estimates), NA_real_, (log_estimates - true_values) / true_values)) |>
  glimpse()

relative_error_CIs <- relative_error |> 
  group_by(parameter, case) |> 
  summarise(median = median(relative_error),
            upper = quantile(relative_error, 0.95, na.rm = T),
            lower = quantile(relative_error, 0.05, na.rm = T)) |> 
  glimpse()


save(relative_error_CIs, file = here::here("data", "HPCC_relative_error_with_CIs_v2.Rdata"))



# code from first HPCC run to calculate prop_converged; not needed --------

# calculate proportion converged ------------------------------------------
# pc = proportion converged
# cases a, b, e, f have 19 params; else 23

param_names <- HPCC_results |> 
  filter(!str_detect(parameter, "cor")) |> 
  filter(!str_detect(parameter, "theta")) |>
  #filter(parameter != "log_sigma_PE") |> # i dont' think the HPCC can recover this small of values
  group_by(parameter) |> 
  tally() |> pull(parameter)

pc_19 <- HPCC_results |> 
  filter(case %in% c("a", "b", "e", "f")) |> 
  mutate(sim = rep(1:4000, each = 19))
pc_23 <- HPCC_results |> 
  filter(!case %in% c("a", "b", "e", "f")) |> 
  mutate(sim = rep(1:4000, each = 23))
pc_full <- bind_rows(pc_19, pc_23)

pc <- pc_full |> 
  filter(parameter %in% param_names) |> 
  group_by(case, sim) |> 
  summarise(nas = if_else(is.na(log_scale_SE), TRUE, FALSE),
            converged = if_else(sum(nas) > 0, FALSE, TRUE)) |> 
  ungroup() |> 
  distinct(case, sim, converged) |> 
  group_by(case) |> 
  summarise(prop_converged = sum(converged)/max(sim)) |> glimpse()

save(pc, file = here::here("data", "HPCC_results_prop_converged.Rdata"))

# old code for calculating prop converged ---------------------------------


# this assumes that if one of the params failed, the others did (i.e. does not assume overlap by taking max)
prop_HPCC_results_converged <- HPCC_results |>
  filter(!str_detect(parameter, "cor")) |>  # these will not have standard errors
  filter(!str_detect(parameter, "theta")) |> 
  mutate(converged = if_else(is.na(log_scale_SE), TRUE, FALSE)) |> # check if a param has no log_scale_SE
  group_by(case, parameter) |> 
  summarise(sum_na_params = sum(converged)/8000) |> # get prop of a param that didn't get a SE
  ungroup() |> 
  group_by(case) |> # now take the highest prop for a param for a case; this is prop that didn't converge
  summarise(prop_converged_by_worst_parameter = (1 - max(sum_na_params)))

prop_HPCC_results_converged |> print(n=Inf)

save_filepath <- here::here("data", "HPCC_results_prop_converged.Rdata")
save(prop_HPCC_results_converged, file = save_filepath)



# scratch -----------------------------------------------------------------

# raw CIs
HPCC_results_CIs <- HPCC_results |> 
  filter(!is.na(log_scale_SE)) |> 
  group_by(parameter, case) |> 
  summarise(median = median(log_estimates),
            lower = quantile(log_estimates, 0.05),
            upper = quantile(log_estimates, 0.95)) |> 
  glimpse()

# compute relative error
HPCC_results_with_true <- HPCC_results |> 
  filter(!parameter %in% c("theta5", "theta6")) |> 
  filter(!str_detect(parameter, "cor")) |>
  left_join(estimates_fit |> 
              mutate(true_estimates = as.numeric(log_estimates)) |> 
              select(true_estimates, parameter), 
            by = "parameter") |> 
  mutate(log_estimates = as.numeric(log_estimates),
         relative_error = (log_estimates - true_estimates) / true_estimates) |> 
  group_by(parameter) |> 
  summarise(median = median(relative_error),
            upper = quantile(relative_error, 0.95),
            lower = quantile(relative_error, 0.05)) |> 
  glimpse()


get_relative_error <- function(results, true_estimates_fit) {
  
  results <- results |> 
    mutate(log_scale_SE = as.numeric(log_scale_SE),
           log_estimates = as.numeric(log_estimates))
  
  cases_not_estimating_thetas <- results |> 
    filter(!case %in% c("c", "d", "g", "h"))
  
  theta5_vec <- tibble(parameter = "theta5", log_estimates = as.character(0))
  theta6_vec <- tibble(parameter = "theta6", log_estimates = as.character(0))
  
  estimates_fit_without_thetas <- estimates_fit |>
    bind_rows(theta5_vec, theta6_vec) |>
    mutate(true_vals = as.numeric(log_estimates)) |>
    dplyr::select(true_vals, parameter)
  
  summary_without_thetas <- cases_not_estimating_thetas |>
    left_join(estimates_fit_without_thetas,
              by = "parameter") |>
    group_by(parameter) |> 
    summarise(relative_error = (log_estimates - true_vals) / true_vals,
              median = median(relative_error),
              upper = quantile(relative_error, 0.95),
              lower = quantile(relative_error, 0.05))
  
  
  if(sum(results$case %in% c("c", "d", "g", "h")) > 1){
    
    cases_estimating_thetas <- results |> 
      filter(case %in% c("c", "d", "g", "h"))
    
    theta5_vec <- tibble(parameter = "theta5", log_estimates = as.character(0.58))
    theta6_vec <- tibble(parameter = "theta6", log_estimates = as.character(0.58))
    
    estimates_fit_with_thetas <- estimates_fit |>
      bind_rows(theta5_vec, theta6_vec) |>
      mutate(true_vals = as.numeric(log_estimates)) |>
      dplyr::select(true_vals, parameter)
    
    summary_with_thetas <- cases_estimating_thetas |> 
      left_join(estimates_fit_with_thetas,
                by = "parameter") |>
      group_by(parameter) |> 
      summarise(relative_error = (log_estimates - true_vals) / true_vals,
                median = median(relative_error),
                upper = quantile(relative_error, 0.95),
                lower = quantile(relative_error, 0.05))
  } 
  
  final_summary <- bind_rows(summary_without_thetas, summary_with_thetas)
  
  return(final_summary)
}

