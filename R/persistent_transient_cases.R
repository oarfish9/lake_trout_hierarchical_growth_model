# Stebbins 1-15-2024
# This script takes the best fitting model and simulates data with different
# levels of persistent and transient error, and then fits the model

library(TMB)
library(tidyverse)

# load helper functions
source(here("helper", "persistent_transient_cases_helper.R"))

# prepare for simulations ---------------------------------------
# load dataset
load(here("data-raw", "laketrout_inc_dat.Rdata"))
# create sim vectors
maxage <- 50
npop_sim <- 6
sim_vectors = sim_block_vectors(maxage, npop_sim)
# load parameter estimates for simulations
load(here("data", "obj_env_last_par_true.Rdata"))

# create simulation object
obj <- create_MVNMVN_obj(idat, sim_vectors)$obj

# cases ------------------------------------------------------------------
case_1 <- run_cases(obj, obj_env_last_par, 1) |> mutate(case = 1)
case_2 <- run_cases(obj, obj_env_last_par, 2) |> mutate(case = 2)
case_3 <- run_cases(obj, obj_env_last_par, 3) |> mutate(case = 3)
case_4 <- run_cases(obj, obj_env_last_par, 4) |> mutate(case = 4)

all_cases <- bind_rows(case_1, case_2, case_3, case_4) |> 
  mutate(case = as.character(case))

save(all_cases, file = here("data", "persistent_transient_sim_results.Rdata"))

