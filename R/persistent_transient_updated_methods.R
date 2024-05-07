# updated persistent and transient simulation methods
# according to feedback from reviewer 4-22-2024

# fit model with only persistent error and save estimates
source(here("R", "fit_model_no_Transient_error.R"))

# fit model with only transient error and save estimates
source(here("R", "fit_model_no_persistent_error.R"))

# now do simulations ------------------------------------------------------

source(here("helper", "persistent_transient_cases_helper.R"))
load(here("data", "obj_env_last_par_no_persistent.Rdata"))
load(here("data", "obj_env_last_par_no_transient.Rdata"))
load(here("data", "obj_env_last_par_true.Rdata"))

# create sim vectors
maxage <- 50
npop_sim <- 6
sim_vectors = sim_block_vectors(maxage, npop_sim)

# create simulation object
obj <- create_MVNMVN_obj(idat, sim_vectors)$obj

# simulate data for all three cases using persistent error only
no_transient_case_2 <- run_cases_updated_methods(obj, obj_env_last_par_no_transient) |> mutate(case = 2)
no_persistent_case_3 <- run_cases_updated_methods(obj, obj_env_last_par_no_persistent) |> mutate(case = 3)


# compare to the simulated data for all three cases using full obj_env_last_par (both error)
load(here("data", "persistent_transient_sim_results.Rdata"))


# plot results ------------------------------------------------------------


# plots
top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")

# updated plot 5-1-2024
correct_sims <- all_cases |> 
  mutate(case = as.numeric(case)) |> 
  filter(case == 1) |> 
  bind_rows(no_transient_case_2) |> 
  bind_rows(no_persistent_case_3) |> 
  filter(age %in% c(4, 15, 40)) |> 
  glimpse()

pers_trans_plot_2024_05_06 <- correct_sims |> 
  mutate(case = as.character(case)) |> 
  mutate(case = case_when(case == "1" ~ "Per + Trans",
                          case == "2" ~ "Per",
                          case == "3" ~ "Trans",
                          TRUE ~ case),
         age_title = paste0("Age = ", age),
         age_title_f = factor(age_title, levels = c("Age = 4", "Age = 15", "Age = 40"))) |> 
  ggplot(aes(x = sim_length, color = case)) + 
  geom_density(size = 1.2) +
  facet_wrap(~age_title_f, nrow = 3, scales = "free") +
  scale_color_manual(values = top_palette) + 
  theme_bw() +
  labs(x = "Simulated True Length",
       y = "Density",
       title = "Density plot of simulated true lengths using estimates from models with different assumptions") +
  labs(color = "Case")

ggsave(pers_trans_plot_2024_05_06, file = here("figures", "pers_trans_plot_2024_05_06.png"))

# compare AIC

load(here("data", "AIC_table.Rdata"))
load(here("data", "AIC_no_persistent_model.Rdata"))
load(here("data", "AIC_no_transient_model.Rdata"))

best_AIC <- AIC_table |> 
  filter(convergence == "YES") |> 
  slice_min(AIC) |> 
  pull(AIC)

AIC_comparison_table_sims <- tibble("model_description" = c("all error", 
                                                            "no transient error",
                                                            "no persistent error"),
                                    "AIC" = c(best_AIC,
                                              AIC_no_transient,
                                              AIC_no_persistent))


# scratch -----------------------------------------------------------------

# from when i was running all cases with obj_env_last_par from a model with no transient,
# obj_env_last_par from a model with no persistent


pers_case_1 <- run_cases_updated_methods(obj, obj_env_last_par_no_transient, 1) |> mutate(case = 1) 
pers_case_2 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 2) |> mutate(case = 2) 
pers_case_3 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 3) |> mutate(case = 3) 
pers_case_4 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 4) |> mutate(case = 4) 

no_transient_sim_cases <- bind_rows(pers_case_1, pers_case_2, pers_case_3, pers_case_4) 
save(no_transient_sim_cases, file = here("data", "persistent_transient_sim_results_no_transient.Rdata"))

# simulate data for all three cases using transient error only
trans_case_1 <- run_cases_updated_methods(obj, obj_env_last_par_no_PE, 1) |> mutate(case = 1) 
trans_case_2 <- run_cases_updated_methods(obj, obj_env_last_par_no_PE, 2) |> mutate(case = 2) 
trans_case_3 <- run_cases_updated_methods(obj, obj_env_last_par_no_PE, 3) |> mutate(case = 3) 
trans_case_4 <- run_cases_updated_methods(obj, obj_env_last_par_no_PE, 4) |> mutate(case = 4) 

no_persistent_sim_cases <- bind_rows(trans_case_1, trans_case_2, trans_case_3, trans_case_4) 
save(no_persistent_sim_cases, file = here("data", "persistent_transient_sim_results_no_persistent.Rdata"))


# bind into one DF with column "source point estimates" 

all_persistent_transient_datasets <- bind_rows(no_persistent_sim_cases |> 
                                                 mutate(estimation_scenario = "no_persistent"),
                                               no_transient_sim_cases |> 
                                                 mutate(estimation_scenario = "no_transient"),
                                               all_cases |> 
                                                 mutate(estimation_scenario = "both",
                                                        case = as.numeric(case)))

# plots
load(here("data", "persistent_transient_sim_results_no_persistent.Rdata"))
load(here("data", "persistent_transient_sim_results_no_transient.Rdata"))
load(here("data", "persistent_transient_sim_results.Rdata"))

no_transient_sim_cases_plot <- no_transient_sim_cases |> 
  mutate(case = as.character(case)) |> 
  filter(age %in% c(4, 15, 40),
         sim_length < 2e05,
         case != "4") |> 
  mutate(case = case_when(case == "1" ~ "Per + Trans",
                          case == "2" ~ "Per",
                          case == "3" ~ "Trans",
                          TRUE ~ case),
         age_title = paste0("Age = ", age),
         age_title_f = factor(age_title, levels = c("Age = 4", "Age = 15", "Age = 40"))) |> 
  ggplot(aes(x = sim_length, color = case)) + 
  geom_density(size = 1.2) +
  facet_wrap(~age_title_f, nrow = 3, scales = "free") +
  scale_color_manual(values = top_palette) + 
  #theme_minimal() +
  theme_bw() +
  labs(x = "Simulated True Length",
       y = "Density",
       title = "Density plot of simulated true lengths using estimates from a model fit with no transient error") +
  labs(color = "Case")

ggsave(no_transient_sim_cases_plot, file = here("figures", "no_transient_sim_cases_plot.png"))

no_persistent_sim_cases_plot <- no_persistent_sim_cases |> 
  mutate(case = as.character(case)) |> 
  filter(age %in% c(4, 15, 40),
         sim_length < 1000,
         case != "4") |> 
  mutate(case = case_when(case == "1" ~ "Per + Trans",
                          case == "2" ~ "Per",
                          case == "3" ~ "Trans",
                          TRUE ~ case),
         age_title = paste0("Age = ", age),
         age_title_f = factor(age_title, levels = c("Age = 4", "Age = 15", "Age = 40"))) |> 
  ggplot(aes(x = sim_length, color = case)) + 
  geom_density(size = 1.2) +
  facet_wrap(~age_title_f, nrow = 3, scales = "free") +
  scale_color_manual(values = top_palette) + 
  #theme_minimal() +
  theme_bw() +
  labs(x = "Simulated True Length",
       y = "Density",
       title = "Density plot of simulated true lengths using estimates from a model fit with no persistent error") +
  labs(color = "Case")

ggsave(no_persistent_sim_cases_plot, file = here("figures", "no_persistent_sim_cases_plot.png"))

