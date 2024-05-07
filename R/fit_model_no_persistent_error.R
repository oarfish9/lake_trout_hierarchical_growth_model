# libraries
library(tidyverse)
library(TMB)
library(here)

# hard code model with no transient error to see if it converges
# May 1, 2024

# helper functions
source(here("helper", "best_fitting_model_helper.R"))

# Compile the model and load in the functions
setwd(here("TMB"))
if(is.loaded("fit_model_no_persistent_error.cpp")){
  dyn.unload("fit_model_no_persistent_error")
}
TMB::compile("fit_model_no_persistent_error.cpp")
dyn.load("fit_model_no_persistent_error")

# load the data
load(here("data-raw", "laketrout_inc_dat.Rdata"))

# Some useful numbers for indexing
nrec <- nrow(idat)

# The data to use in the model
data <- list(
  nrec = nrec, 
  L = idat$L,
  age = idat$Age
)

# List the parameters and give starting values
parameters <- list(
  log_sigma = log(40), # log(40) worked, so did a log(10) and log(500) combo
  log_Linf_hyper = log(300), # 300
  log_K_hyper = log(0.1),
  log_L1_hyper = log(120), # 120
  log_sigma_PE = log(1), # log(1) worked
  log_PE = rep(0, nrec) # vector for storing process error, no space for cases where age=0
)

map_par <- list()


# Define the random effects
RE <- c('log_PE')

# Generate the objective function and gradient
obj <- MakeADFun(
  data = data,
  parameters = parameters,
  DLL = 'fit_model_no_persistent_error',
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
AIC_no_persistent <- (2*p) - (2*-NLL) # switch sign of NLL, AIC calculation

# for simulations
obj_env_last_par_no_persistent <- obj$env$last.par

save(obj_env_last_par_no_persistent, file = here("data", "obj_env_last_par_no_persistent.Rdata"))
save(AIC_no_persistent, file = here("data", "AIC_no_persistent_model.Rdata"))

