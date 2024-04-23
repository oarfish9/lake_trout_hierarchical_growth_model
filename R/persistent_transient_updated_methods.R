# updated persistent and transient simulation methods
# according to feedback from reviewer 4-22-2024

# fit model with only persistent error and save estimates

# fit model with only transient error and save estimates

# now simulate data using persistent point estimates

# now simulate data with no transient error using transient point estimates


# libraries
library(tidyverse)
library(TMB)
library(here)

# helper functions
source(here("helper", "best_fitting_model_helper.R"))

# Compile the model and load in the functions
setwd(here("TMB"))
if(is.loaded("ltSup_MVN_BVN_v2.cpp")){
  dyn.unload("ltSup_MVN_BVN_v2.cpp")
}
TMB::compile("ltSup_MVN_BVN_v2.cpp")
dyn.load("ltSup_MVN_BVN_v2")

# load the data
load(here("data-raw", "laketrout_inc_dat.Rdata"))

fit_model_per_or_trans <- function(per_or_trans) {
  # Some useful numbers for indexing
  nrec <- nrow(idat)
  nind <- length(unique(idat$idx))
  npop <- length(unique(idat$pop_idx))
  pop_idx_ind <- idat[!duplicated(idat$idx),]$pop_idx #length of nind, values of npop for indexing
  
  # create sim vectors
  maxage <- 50
  npop_sim <- 6
  sim_vectors = sim_block_vectors(maxage, npop_sim)
  
  # The data to use in the model
  data <- list(
    ind_idx = idat$idx,
    pop_idx_ind = pop_idx_ind,
    ind_idx_sim = sim_vectors$ind_idx_sim,
    pop_idx_ind_sim = sim_vectors$pop_idx_ind_sim,
    npop = npop,
    npop_sim = npop_sim,
    nrec = nrec, 
    nind = nind,
    L = idat$L,
    age = idat$Age,
    maxage = maxage,
    age_sim = sim_vectors$age
  )
  
  # List the parameters and give starting values
  parameters <- list(
    log_sigma = log(1),
    log_sigma_log_Linf_devs = log(30),
    log_sigma_log_K_devs = log(0.05), # log of sd of log K deviations
    log_sigma_log_L1_devs = log(0.05),# need a starting value
    log_sigma_log_Linf_hyper_devs = log(30),
    log_sigma_log_K_hyper_devs = log(0.05),
    log_sigma_log_L1_hyper_devs = log(0.05),
    log_Linf_hyper = log(644),
    log_K_hyper = log(0.10),
    log_L1_hyper = log(100),
    theta1 = 0, # set to 0, = no priors; theta1 = Linf/K
    theta2 = 0,
    theta3 = 0,
    theta4 = 0,
    log_sigma_PE = log(4), # sigma for process error estimation
    log_Linf_devs = rep(0, nind), # log of deviations for each individual fish
    log_K_devs = rep(0, nind),
    log_L1_devs = rep(0, nind),
    log_Linf_hyper_devs = rep(0, npop), # log of deviations for population dist.
    log_K_hyper_devs = rep(0, npop),
    log_L1_hyper_devs = rep(0, npop),
    log_PE = rep(0, nrec) # vector for storing process error, no space for cases where age=0
  )
  
  # no need to map pars
  if(per_or_trans == "per") {
    # do not estimate transient error
    map_par <- list(log_sigma_PE = factor(NA))
  } else {
    # do not estimate persistent error
    map_par <- list(log_sigma_log_Linf_devs = factor(NA),
                    log_sigma_log_K_devs =factor(NA),
                    log_sigma_log_L1_devs = factor(NA),
                    log_sigma_log_Linf_hyper_devs = factor(NA),
                    log_sigma_log_L1_hyper_devs = factor(NA),
                    log_sigma_log_K_hyper_devs = factor(NA),
                    theta1 = factor(NA),
                    theta2 = factor(NA),
                    theta3 = factor(NA),
                    theta4 = factor(NA))
  }
  
  # Define the random effects
  RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
          'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')
  
  # Generate the objective function and gradient
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = 'ltSup_MVN_BVN_v2',
    random = RE,
    map = map_par
  )
  
  # run the model
  opt <- nlminb(start=obj$par, objective=obj$fn, 
                gradient=obj$gr)
  
  # report out
  rep <- obj$report()
  sdr <- sdreport(obj)
  
  # format output into tibble
  estimates_fit <- format_output(sdr, rep, ind_level = "MVN",
                                 pop_level = "BVN") |> 
    print(n = Inf)  
  
  NLL <- rep$f
  p <- length(parameters) - length(map_par) # number of parameters
  AIC <- (2*p) - (2*-NLL) # switch sign of NLL, AIC calculation
  
  # for simulations
  obj_env_last_par <- obj$env$last.par
  
  
  return(list("rep" = rep,
              "sdr" = sdr,
              "estimates_fit" = estimates_fit,
              "AIC" = AIC,
              "obj_env_last_par" = obj_env_last_par))
}

point_estimates_no_transient <- fit_model_per_or_trans("per")
point_estimates_no_persistent <- fit_model_per_or_trans("trans")

save(point_estimates_no_transient, file = here("data", "point_estimates_for_per_trans_no_transient.Rdata"))
save(point_estimates_no_persistent, file = here("data", "point_estimates_for_per_trans_no_persistent.Rdata"))


# now do simulations ------------------------------------------------------

source(here("helper", "persistent_transient_cases_helper.R"))

# create sim vectors
maxage <- 50
npop_sim <- 6
sim_vectors = sim_block_vectors(maxage, npop_sim)

# create simulation object
obj <- create_MVNMVN_obj(idat, sim_vectors)$obj

# simulate data for all three cases using persistent error only
pers_case_1 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 1) |> mutate(case = 1) 
pers_case_2 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 2) |> mutate(case = 1) 
pers_case_3 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 3) |> mutate(case = 1) 
pers_case_4 <- run_cases_updated_methods(obj, point_estimates_no_transient$obj_env_last_par, 4) |> mutate(case = 1) 

no_transient_sim_cases <- bind_rows(pers_case_1, pers_case_2, pers_case_3, pers_case_4) 
save(no_transient_sim_cases, file = here("data", "persistent_transient_sim_results_no_transient.Rdata"))

# simulate data for all three cases using transient error only
trans_case_1 <- run_cases_updated_methods(obj, point_estimates_no_persistent$obj_env_last_par, 1) |> mutate(case = 1) 
trans_case_2 <- run_cases_updated_methods(obj, point_estimates_no_persistent$obj_env_last_par, 2) |> mutate(case = 1) 
trans_case_3 <- run_cases_updated_methods(obj, point_estimates_no_persistent$obj_env_last_par, 3) |> mutate(case = 1) 
trans_case_4 <- run_cases_updated_methods(obj, point_estimates_no_persistent$obj_env_last_par, 4) |> mutate(case = 1) 

no_persistent_sim_cases <- bind_rows(trans_case_1, trans_case_2, trans_case_3, trans_case_4) 
save(no_persistent_sim_cases, file = here("data", "persistent_transient_sim_results_no_persistent.Rdata"))

# compare to the simulated data for all three cases using full obj_env_last_par (both error)
load(here("data", "persistent_transient_sim_results.Rdata"))

# compare AIC


