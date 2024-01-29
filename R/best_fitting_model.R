# Stebbins 1-15-2024
# This takes the correlation structure identified by evaluate_alternate_models.R
# and fits the model, saving the output

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
map_par <- list()

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
estimates_fit <- format_output(sdr, rep, MVN_or_BVN = "BVN") |> 
  print(n = Inf)

write.csv(estimates_fit, here("data", "true_estimates.csv"))
obj_env_last_par <- obj$env$last.par
save(estimates_fit, file = "true_estimates.Rdata")
save(obj_env_last_par, file = here("data", "obj_env_last_par_true.Rdata"))

# save parameters for use in simulation_code_for_HPCC.R
# insert placeholders for theta5 and theta6
original_pars_1 <- obj_env_last_par[1:14]
original_pars_placeholder <- c("theta5" = 0,
                               "theta6" = 0)
original_pars_2 <- obj_env_last_par[15:length(obj_env_last_par)]
original_pars_for_sim <- c(original_pars_1, original_pars_placeholder,
                           original_pars_2)
save(original_pars_for_sim, file = here("data", "original_pars_for_sim.Rdata"))

