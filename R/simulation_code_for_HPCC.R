# load libraries
rm(list=ls())
library(tidyverse)
library(TMB)
library(ggplot2)
library(lattice)
library(here)

# functions ---------------------------------------------------------------
source(here("helper", "simulation_code_for_HPCC_helper.R"))

# create vectors for simulation block -------------------------------------
npop_sim = 6
maxage = 50
sim_vectors <- sim_block_vectors(maxage = 50, npop_sim = 6)

load(here("data", "original_pars_for_sim.Rdata"))
# prep to run simulations with 6 populations -------------------------------------------------

nsims <- 20
npop_sim <- 6
maxage <- 50 # maximum age in years
sim_vectors <- sim_block_vectors(maxage, npop_sim) # create sim vectors
obj_sim <- create_MVNMVN_obj(idat, sim_vectors)$obj # need obj that will estimate theta5 and theta6


# case A  ----------------------------------------------------------------

results_case_a <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0, theta6 = 0, fixed = TRUE)

case_a_output_clean <- results_case_a$raw_output


# run case B -------------------------------------------------------------

results_case_b <- run_simulations(nsims = nsims, sim_vectors = sim_vectors, 
                                  maxage = maxage, obj_sim = obj_sim,
                                  npop_sim = npop_sim, 
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0.58, theta6 = 0.58, fixed = TRUE)

case_b_output_clean <- results_case_b$raw_output




# run case C --------------------------------------------------------------
results_case_c <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors,
                                  maxage = maxage,
                                  obj_sim = obj_sim,
                                  npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0, theta6 = 0, 
                                  fixed = FALSE)

case_c_output_clean <- results_case_c$raw_output                               

# run case D --------------------------------------------------------------

results_case_d <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0.58, theta6 = 0.58, fixed = FALSE)

case_d_output_clean <- results_case_d$raw_output


# prep to run simulations with 50 populations -------------------------------------------------  
npop_sim = 50
maxage = 50 # maximum age in years
sim_vectors <- sim_block_vectors(maxage, npop_sim) # need to re run for some reason
obj_sim <- create_MVNMVN_obj(idat, sim_vectors)$obj

# run case E --------------------------------------------------------------
results_case_e <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0, theta6 = 0, fixed = TRUE)

case_e_output_clean <- results_case_e$raw_output

# run case F --------------------------------------------------------------

results_case_f <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0.58, theta6 = 0.58, fixed = TRUE)

case_f_output_clean <- results_case_f$raw_output


# run case G --------------------------------------------------------------

results_case_g <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0, theta6 = 0, fixed = FALSE)

case_g_output_clean <- results_case_g$raw_output

# run case H --------------------------------------------------------------

results_case_h <- run_simulations(nsims = nsims,
                                  sim_vectors = sim_vectors, maxage = maxage,
                                  obj_sim = obj_sim, npop_sim = npop_sim,
                                  original_parameter_estimates = original_pars_for_sim,
                                  theta5 = 0.58, theta6 = 0.58, fixed = FALSE)

case_h_output_clean <- results_case_h$raw_output
