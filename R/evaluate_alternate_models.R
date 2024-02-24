# Stebbins
# 1-15-2024
# This script evaluates versions of the hierarchical von Bertalanffy model
# with different correlation structures for the growth parameters L infinity,
# K, and L1 at the individual and population levels.

# load libraries
library(tidyverse)
library(TMB)
library(ggplot2)
library(lattice)
library(here)

# source helper functions
source(here("helper", "evaluate_alternate_models_helper.R"))

# create vectors for simulation block -------------------------------------
npop_sim = 6 # six simulated populations
maxage = 50 # max age of fish in simulated populations
sim_vectors <- sim_block_vectors(maxage = 50, npop_sim = 6) # create vectors for the simulation block

# function to fit and compare growth models -------------------------------

# pass in what you want the distribution to be at the individual level (multivariate normal, bivariate normal (Linf and K), or univariate (no correlations))
evaluate_alternate_models <- function(ind_level = c("MVN", "BVN", "UVN"),
                                      pop_level = c("MVN", "BVN", "UVN")) {
  
  # compile the model and load in the functions
  setwd(here("TMB"))
  if(is.loaded("ltSup_MVN_MVN_setnpopsim.cpp")){
    dyn.unload("ltSup_MVN_MVN_setnpopsim.cpp")
  }
  TMB::compile("ltSup_MVN_MVN_setnpopsim.cpp")
  dyn.load("ltSup_MVN_MVN_setnpopsim")
  
  # load the data
  load(here("data-raw", "laketrout_inc_dat.Rdata"))
  
  # Some useful numbers for indexing
  nrec <- nrow(idat)
  nind <- length(unique(idat$idx))
  npop <- length(unique(idat$pop_idx))
  pop_idx_ind <- idat[!duplicated(idat$idx),]$pop_idx #length of nind, values of npop for indexing
  
  # the data to use in the model
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
  
  # list the parameters and give starting values
  parameters <- list(
    log_sigma = log(1),
    log_sigma_log_Linf_devs = log(30),
    log_sigma_log_K_devs = log(0.05), # log of sd of log K deviations
    log_sigma_log_L1_devs = log(0.05), # uninformative
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
    theta5 = 0,
    theta6 = 0,
    log_sigma_PE = log(4), # sigma for process error estimation
    log_Linf_devs = rep(0, nind), # log of deviations for each individual fish
    log_K_devs = rep(0, nind),
    log_L1_devs = rep(0, nind),
    log_Linf_hyper_devs = rep(0, npop), # log of deviations for population dist.
    log_K_hyper_devs = rep(0, npop),
    log_L1_hyper_devs = rep(0, npop),
    log_PE = rep(0, nrec) # vector for storing process error, no space for cases where age=0
  )
  
  # map correlation parameters according to arguments
  # parameters in map_par will be estimated
  if(ind_level == "MVN" & pop_level == "MVN") {
    map_par <- list()
  } else if(ind_level == "MVN" & pop_level == "BVN"){
    map_par <- list(
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  } else if(ind_level == "MVN" & pop_level == "UVN"){
    map_par <- list(
      theta4 = factor(NA), # Linf/K pop
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  } else if(ind_level == "BVN" & pop_level == "MVN"){
    map_par <- list(
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA)  # K/L1 ind
    )
  } else if(ind_level == "BVN" & pop_level == "BVN"){
    map_par <- list(
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA), # K/L1 ind
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  } else if(ind_level == "BVN" & pop_level == "UVN"){
    map_par <- list(
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA), # K/L1 ind
      theta4 = factor(NA), # Linf/K pop
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  } else if(ind_level == "UVN" & pop_level == "MVN"){
    map_par <- list(
      theta1 = factor(NA), # Linf/K ind
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA)  # K/L1 ind
    )
  } else if(ind_level == "UVN" & pop_level == "BVN"){
    map_par <- list(
      theta1 = factor(NA), # Linf/K ind
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA), # K/L1 ind
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  } else if(ind_level == "UVN" & pop_level == "UVN"){
    map_par <- list(
      theta1 = factor(NA), # Linf/K ind
      theta2 = factor(NA), # Linf/L1 ind
      theta3 = factor(NA), # K/L1 ind
      theta4 = factor(NA), # Linf/K pop
      theta5 = factor(NA), # Linf/L1 pop
      theta6 = factor(NA)  # K/L1 pop
    )
  }
  
  # random effects (individua- and population-level deviations from the mean value)
  RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
          'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')
  
  # generate the objective function and gradient
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = "ltSup_MVN_MVN_setnpopsim",
    random = RE,
    map = map_par
  )
  
  # run the model
  opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr)
  
  # report out
  rep <- obj$report()
  sdr <- sdreport(obj)
  
  # format output into tibble
  estimates_fit <- format_output(sdr, rep, 
                                 ind_level = ind_level, 
                                 pop_level = pop_level) |> 
    mutate(log_scale_SE = if_else(log_scale_SE == "NaN", NA_character_, log_scale_SE),
           log_scale_SE = as.numeric(log_scale_SE),
           log_estimates = round(as.numeric(log_estimates), 5))
  
  # calculate AIC
  NLL <- rep$f
  p <- length(parameters) - length(map_par) # number of parameters
  AIC <- (2*p) - (2*-NLL) # switch sign of NLL, AIC calculation
  
  return(list("estimates_fit" = estimates_fit,
              "AIC" = AIC))
  
}

# describe models ---------------------------------------------------------

model_descriptions <- tibble(model = 1:9,
                             letter = LETTERS[1:9],
                             ind_dist = c(rep("MVN", 3), 
                                          rep("BVN", 3), 
                                          rep("UVN", 3)),
                             pop_dist = c(rep(c("MVN", "BVN", "UVN"), 3)))

save(model_descriptions, file = here("data", "model_descriptions_table.Rdata"))


# run the models ----------------------------------------------------------

model_1 <- evaluate_alternate_models(ind_level = "MVN", pop_level = "MVN") |> 
  get_estimates() |> 
  mutate(model = 1)

model_2 <- evaluate_alternate_models(ind_level = "MVN", pop_level = "BVN") |> 
  get_estimates() |> 
  mutate(model = 2)

model_3 <- evaluate_alternate_models(ind_level = "MVN", pop_level = "UVN") |> 
  get_estimates() |> 
  mutate(model = 3)

model_4 <- evaluate_alternate_models(ind_level = "BVN", pop_level = "MVN") |> 
  get_estimates() |> 
  mutate(model = 4)

model_5 <- evaluate_alternate_models(ind_level = "BVN", pop_level = "BVN") |> 
  get_estimates() |> 
  mutate(model = 5)

model_6 <- evaluate_alternate_models(ind_level = "BVN", pop_level = "UVN") |> 
  get_estimates() |> 
  mutate(model = 6)

model_7 <- evaluate_alternate_models(ind_level = "UVN", pop_level = "MVN") |> 
  get_estimates() |> 
  mutate(model = 7)

model_8 <- evaluate_alternate_models(ind_level = "UVN", pop_level = "BVN") |> 
  get_estimates() |> 
  mutate(model = 8)

model_9 <- evaluate_alternate_models(ind_level = "UVN", pop_level = "UVN") |> 
  get_estimates() |> 
  mutate(model = 9)

all_results <- bind_rows(model_1, model_2, model_3, model_4, model_5,
                         model_6, model_7, model_8, model_9)

save(all_results, file = here("data", "model_evaluation_results.Rdata"))


# check for convergence
non_converged_models <- all_results |> 
  filter(!str_detect(parameter, "cor")) |> 
  filter(is.na(log_scale_SE)) |> 
  distinct(model) |> 
  pull(model)

# gather AIC
AIC_table <- all_results |> 
  distinct(model, AIC) |> 
  mutate(convergence = case_when(model %in% non_converged_models ~ "NO",
                                 TRUE ~ "YES"),
         letter = LETTERS[model])

save(AIC_table, file = here::here("data", "AIC_table.Rdata"))

# calculate probability intervals ------------------------------------------

probability_intervals <- all_results |> 
  mutate(letter = LETTERS[model]) |> 
  filter(letter == "B") |> 
  mutate(log_estimates = round(log_estimates, 4),
         estimates = round(estimates, 4),
         CI_LB = round(exp(log_estimates - (2 * log_scale_SE)), 3),
         CI_UB = round(exp(log_estimates + (2 * log_scale_SE)), 3)) |> 
  select(parameter,
         log_estimates,
         estimates,
         log_scale_SE,
         CI_LB, CI_UB)
save(probability_intervals, file = here::here("data", "best_fitting_model_with_PIs.Rdata"))

