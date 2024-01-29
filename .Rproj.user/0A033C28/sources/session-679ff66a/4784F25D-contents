# Stebbins 1-15-2024 
# This script runs the supplementary models: biphasic and sex-specific
library(tidyverse)
library(TMB)

# functions ---------------------------------------------------------------
source(here("helper", "biphasic_sexspecific_helper.R"))

# run model and get obj_env_last_par --------------------------------------
maxage <- 50
npop_sim <- 6
sim_vectors_biphasic <- sim_block_vectors(maxage, npop_sim)
obj_biphasic <- create_biphasic_obj(sim_vectors_biphasic)$obj

opt <- nlminb(start=obj_biphasic$par, objective=obj_biphasic$fn, 
              gradient=obj_biphasic$gr)

rep <- obj_biphasic$report()
sdr <- sdreport(obj_biphasic)

biphasic_estimates <- format_output_biphasic(sdr, rep)
biphasic_obj_env_last_par <- obj_biphasic$env$last.par

save(biphasic_estimates, file = here("data", "biphasic_estimates.Rdata"))
save(biphasic_obj_env_last_par, file = here("data", "biphasic_obj_env_last_par.Rdata"))

# now run simulations (biphasic) -----------------------------------------------------

# cases ------------------------------------------------------------------
case_1 <- run_biphasic_cases(obj_biphasic, biphasic_obj_env_last_par, 1) |> mutate(case = 1)
case_2 <- run_biphasic_cases(obj_biphasic, biphasic_obj_env_last_par, 2) |> mutate(case = 2)
case_3 <- run_biphasic_cases(obj_biphasic, biphasic_obj_env_last_par, 3) |> mutate(case = 3)
case_4 <- run_biphasic_cases(obj_biphasic, biphasic_obj_env_last_par, 4) |> mutate(case = 4)

all_cases_biphasic <- bind_rows(case_1, case_2, case_3, case_4) |> 
  mutate(case = as.character(case))

save(all_cases_biphasic, file = here("data", "persistent_transient_sim_results_biphasic.Rdata"))

# run model (sex specific) ------------------------------------------------
maxage = 50
npop_sim = 6
sim_vectors_sex_specific <- sim_block_vectors(maxage = 50, npop_sim = 6)

# fit MVN_BVN model to original data -----------------------------------------------------------
load(here("data-raw", "incdat_2022_biodat.Rdata"))
male_dat <- biodat_save |> 
  mutate(idx2 = cumsum(!duplicated(idx))-1,
         idx3 = idx2 + 1)

obj_sex_specific <- create_BVN_MVN_obj(sim_vectors_sex_specific, male_dat)$obj
opt_sex_specific <- nlminb(start=obj_sex_specific$par, objective=obj_sex_specific$fn, 
                           gradient=obj_sex_specific$gr)
rep_sex_specific <- obj_sex_specific$report()
sdr_sex_specific <- sdreport(obj_sex_specific)

sex_specific_estimates <- format_output(sdr_sex_specific, rep_sex_specific, "BVN")
sex_specific_obj_env_last_par <- obj_sex_specific$env$last.par

save(sex_specific_estimates, file = here("data", "sex_specific_estimates.Rdata"))
save(sex_specific_obj_env_last_par, file = here("data", "sex_specific_obj_env_last_par.Rdata"))

# run simulations (sex specific) ------------------------------------------
case_1 <- run_cases(obj_sex_specific, sex_specific_obj_env_last_par, 1) |> mutate(case = 1)
case_2 <- run_cases(obj_sex_specific, sex_specific_obj_env_last_par, 2) |> mutate(case = 2)
case_3 <- run_cases(obj_sex_specific, sex_specific_obj_env_last_par, 3) |> mutate(case = 3)
case_4 <- run_cases(obj_sex_specific, sex_specific_obj_env_last_par, 4) |> mutate(case = 4)

all_cases_sex_specific <- bind_rows(case_1, case_2, case_3, case_4) |> 
  mutate(case = as.character(case))

save(all_cases_sex_specific, file = here("data", "persistent_transient_sim_results_sex_specific.Rdata"))
