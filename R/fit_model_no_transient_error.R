# libraries
library(tidyverse)
library(TMB)
library(here)

# hard code model with no persistent error to see if it converges
# May 6, 2024

# helper functions
source(here("helper", "best_fitting_model_helper.R"))

# Compile the model and load in the functions
setwd(here("TMB"))
if(is.loaded("fit_model_no_transient.cpp")){
  dyn.unload("fit_model_no_transient")
}
TMB::compile("fit_model_no_transient.cpp")
dyn.load("fit_model_no_transient")

# load the data
load(here("data-raw", "laketrout_inc_dat.Rdata"))

# Some useful numbers for indexing
nrec <- nrow(idat)
nind <- length(unique(idat$idx))
npop <- length(unique(idat$pop_idx))
pop_idx_ind <- idat[!duplicated(idat$idx),]$pop_idx #length of nind, values of npop for indexing

# The data to use in the model
data <- list(
  ind_idx = idat$idx,
  pop_idx_ind = pop_idx_ind,
  npop = npop,
  nrec = nrec, 
  nind = nind,
  L = idat$L,
  age = idat$Age
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
  log_Linf_devs = rep(0, nind), # log of deviations for each individual fish
  log_K_devs = rep(0, nind),
  log_L1_devs = rep(0, nind),
  log_Linf_hyper_devs = rep(0, npop), # log of deviations for population dist.
  log_K_hyper_devs = rep(0, npop),
  log_L1_hyper_devs = rep(0, npop)
)

map_par <- list()


# Define the random effects
RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
        'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs')

# Generate the objective function and gradient
obj <- MakeADFun(
  data = data,
  parameters = parameters,
  DLL = 'fit_model_no_transient',
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
AIC_no_transient <- (2*p) - (2*-NLL) # switch sign of NLL, AIC calculation

# for simulations
obj_env_last_par_no_transient <- obj$env$last.par

save(obj_env_last_par_no_transient, file = here("data", "obj_env_last_par_no_transient.Rdata"))
save(AIC_no_transient, file = here("data", "AIC_no_transient_model.Rdata"))

