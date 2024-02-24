# run best fitting model with and without persistent and transient error
# to add to AIC table

# this is not correct - do not use
# setup is the same -------------------------------------------------------

# libraries
library(tidyverse)
library(TMB)

# helper functions
source(here("helper", "best_fitting_model_helper.R"))

# Compile the model and load in the functions
setwd(here("TMB"))
if(is.loaded("ltSup_MVN_BVN_v2")){
  dyn.unload("ltSup_MVN_BVN_v2")
}
TMB::compile("ltSup_MVN_BVN_v2.cpp")
dyn.load("ltSup_MVN_BVN_v2")

# load the data
load(here("data-raw", "laketrout_inc_dat.Rdata"))

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


# fit with no persistent error ---------------------------------------------
# List the parameters and give starting values
parameters_transient <- list(
  log_sigma = log(1),
  log_sigma_log_Linf_devs = 0,
  log_sigma_log_K_devs = 0, # log of sd of log K deviations
  log_sigma_log_L1_devs = 0,# need a starting value
  log_sigma_log_Linf_hyper_devs = 0,
  log_sigma_log_K_hyper_devs = 0,
  log_sigma_log_L1_hyper_devs = 0,
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

map_par_transient <- list(
  log_sigma_log_Linf_devs = factor(NA),
  log_sigma_log_K_devs = factor(NA),
  log_sigma_log_L1_devs = factor(NA),
  log_sigma_log_Linf_hyper_devs = factor(NA),
  log_sigma_log_K_hyper_devs = factor(NA),
  log_sigma_log_L1_hyper_devs = factor(NA),
  theta1 = factor(NA),
  theta2 = factor(NA),
  theta3 = factor(NA),
  theta4 = factor(NA)
)

# list random effects
RE <- c("log_PE")


# Generate the objective function and gradient
obj <- MakeADFun(
  data = data,
  parameters = parameters_transient,
  DLL = 'ltSup_MVN_BVN_v2',
  random = RE,
  map = map_par_transient
)

# run the model
opt <- nlminb(start=obj$par, objective=obj$fn, 
              gradient=obj$gr)

# report out
rep_transient <- obj$report()
sdr_transient <- sdreport(obj)

# fit with no transient error ---------------------------------------------
map_par_persistent <- list(
  log_sigma_PE = factor(NA)
)

# list random effects
RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
        'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')


# Generate the objective function and gradient
obj <- MakeADFun(
  data = data,
  parameters = parameters,
  DLL = 'ltSup_MVN_BVN_v2',
  random = RE,
  map = map_par_persistent
)

# run the model
opt <- nlminb(start=obj$par, objective=obj$fn, 
              gradient=obj$gr)

# report out
rep_persistent <- obj$report()
sdr_persistent <- sdreport(obj)


# calculate AIC -----------------------------------------------------------

NLL_transient <- rep_transient$f
p_transient <- length(parameters) - length(map_par_transient) # number of parameters
AIC_transient <- (2*p_transient) - (2*-NLL_transient) # switch sign of NLL, AIC calculation


NLL_persistent <- rep_persistent$f
p_persistent <- length(parameters) - length(map_par_persistent) # number of parameters
AIC_persistent <- (2*p_persistent) - (2*-NLL_persistent) # switch sign of NLL, AIC calculation
