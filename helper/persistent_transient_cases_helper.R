# Stebbins 1-15-2024
# this contains helper functions for running the script 
# persistent_transient_cases.R

# create age, ind_idx, pop_idx, and pop_ind_idx vectors for simulation block
sim_block_vectors <- function(maxage, npop_sim){
  
  # initialize
  nrec <- maxage * ((maxage + 1)/2)
  nrec_sim <- (maxage * ((maxage + 1)/2)) * npop_sim
  
  ages <- list()
  ind_idx_sim_1 <- list()
  pop_idx_sim <- list()
  
  # age and ind_idx_sim vectors
  for(i in 1:maxage){
    ages[[i]] <- seq(1:i)
    ind_idx_sim_1[[i]] <- rep(i, i)
  }
  ages <- unlist(ages)
  ind_idx_sim_1 <- unlist(ind_idx_sim_1)
  
  age <- rep(ages, npop_sim)
  
  tmp <- list()
  tmp[[1]] <- ind_idx_sim_1
  for(i in 2:npop_sim){
    tmp[[i]] <- rep((i-1) * maxage + ind_idx_sim_1)
  }
  ind_idx_sim <- unlist(tmp)-1
  nind_sim = length(unique(ind_idx_sim))
  
  # now create pop_idx_sim 
  for(i in 1:npop_sim){
    pop_idx_sim[[i]] <- rep(i, nrec)
  }
  pop_idx_sim <- unlist(pop_idx_sim)
  
  sim_age <- data.frame(age, ind_idx_sim, pop_idx_sim)
  
  pop_idx_ind_sim <- unique(sim_age[c('ind_idx_sim', 'pop_idx_sim')])$pop_idx_sim - 1
  
  return(list("age" = age,
              "ind_idx_sim" = ind_idx_sim,
              "pop_idx_sim" = pop_idx_sim,
              "pop_idx_ind_sim" = pop_idx_ind_sim,
              "nrec_sim" = nrec_sim,
              "nind_sim" = nind_sim))
  
}

# function for creating obj$sim for MVN_MVN
create_MVNMVN_obj <- function(idat, sim_vectors) {
  
  # Compile the model and load in the functions
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
  
  # uncomment thetas for a univariate distribution
  map_par <- list()
  
  # Define the random effects
  RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
          'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')
  
  # Generate the objective function and gradient
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = "ltSup_MVN_MVN_setnpopsim",
    random = RE,
    map = map_par
  )
  
  return(list("obj" = obj))
}

# nested function
simulate_data <- function(obj, parameter_list) {
  
  sim_df <- tibble("age" = 0,
                   "sim_length" = 0)
  
  for(i in 1:1000){
    print(paste0("iter: ", i))
    simulate <- obj$simulate(parameter_list, complete = T)
    temp_df <- tibble("age" = simulate$age_sim,
                      "sim_length" = simulate$sim_L_hat,
                      "iter" = i)
    sim_df <- bind_rows(sim_df, temp_df)
  }
  sim_df <- sim_df[-1,]
  return(sim_df)
}

# function for persistent vs transient simulations
run_cases <- function(obj, obj_env_last_par, case) {
  
  # insert placeholders for theta5 and theta6
  original_pars_1 <- obj_env_last_par[1:14]
  original_pars_placeholder <- c("theta5" = 0,
                                 "theta6" = 0)
  original_pars_2 <- obj_env_last_par[15:length(obj_env_last_par)]
  obj_env_last_par <- c(original_pars_1, original_pars_placeholder, 
                        original_pars_2)
  
  if(case == 1) {
    output_data <- simulate_data(obj, obj_env_last_par)
    
  } else if(case == 2) {
    par_list_case_2 <- obj_env_last_par
    par_list_case_2["log_sigma_PE"] <- log(0) # no process error
    output_data <- simulate_data(obj, par_list_case_2)
    
  } else if(case == 3) {
    par_list_case_3 <- obj_env_last_par
    par_list_case_3["log_sigma_log_Linf_devs"] <- log(0) # no persistent error
    par_list_case_3["log_sigma_log_K_devs"] <- log(0)
    par_list_case_3["log_sigma_log_L1_devs"] <- log(0)
    par_list_case_3["log_sigma_log_Linf_hyper_devs"] <- log(0)
    par_list_case_3["log_sigma_log_K_hyper_devs"] <- log(0)
    par_list_case_3["log_sigma_log_L1_hyper_devs"] <- log(0)
    
    output_data <- simulate_data(obj, par_list_case_3)
  } else if(case == 4) {
    par_list_case_4 <- obj_env_last_par
    par_list_case_4["log_sigma_PE"] <- log(0) # no process error
    par_list_case_4["log_sigma_log_Linf_devs"] <- log(0) # no persistent error
    par_list_case_4["log_sigma_log_K_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_devs"] <- log(0)
    par_list_case_4["log_sigma_log_Linf_hyper_devs"] <- log(0)
    par_list_case_4["log_sigma_log_K_hyper_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_hyper_devs"] <- log(0)
    
    output_data <- simulate_data(obj, par_list_case_4)
    
  } else {
    error("provide a value for case that is either 1, 2, 3, or 4")
  }
  return(output_data)
}


# function for persistent vs transient simulations (updated methods)
run_cases_updated_methods <- function(obj, obj_env_last_par, case) {
  
  if(length(obj_env_last_par) == 6479) {
    # fill in persistent error
    original_pars <- obj_env_last_par[1]
    fill_in_pers_parameters <- c("log_sigma_log_Linf_devs" = 0,
                                 "log_sigma_log_K_devs" = 0,
                                 "log_sigma_log_L1_devs" = 0,
                                 "log_sigma_log_Linf_hyper_devs" = 0,
                                 "log_sigma_log_K_hyper_devs" = 0,
                                 "log_sigma_log_L1_hyper_devs" = 0)
    original_pars_pt_2 <- c(obj_env_last_par[2:4])
    fill_in_pers_parameters_2 <- c("theta1" = 0,
                                   "theta2" = 0,
                                   "theta3" = 0,
                                   "theta4" = 0,
                                   "theta5" = 0,
                                   "theta6" = 0)
    original_pars_end <- obj_env_last_par[5:length(obj_env_last_par)]
    obj_env_last_par <- c(original_pars, fill_in_pers_parameters, original_pars_pt_2,
                          fill_in_pers_parameters_2, original_pars_end)
  } 
  
  if(length(obj_env_last_par) == 6488) {
    # fill in transient error and thetas
    original_pars <- obj_env_last_par[1:14]
    fill_in_trans_parameters <- c("theta5" = 0,
                                  "theta6" = 0,
                                  "log_sigma_PE" = 0)
    original_pars_end <- obj_env_last_par[15:length(obj_env_last_par)]
    obj_env_last_par <- c(original_pars, fill_in_trans_parameters, original_pars_end)
  }

  if(case == 1) {
    output_data <- simulate_data(obj, obj_env_last_par)

  } else if(case == 2) {
    par_list_case_2 <- obj_env_last_par
    par_list_case_2["log_sigma_PE"] <- log(0) # no process error
    output_data <- simulate_data(obj, par_list_case_2)

  } else if(case == 3) {
    par_list_case_3 <- obj_env_last_par
    par_list_case_3["log_sigma_log_Linf_devs"] <- log(0) # no persistent error
    par_list_case_3["log_sigma_log_K_devs"] <- log(0)
    par_list_case_3["log_sigma_log_L1_devs"] <- log(0)
    par_list_case_3["log_sigma_log_Linf_hyper_devs"] <- log(0)
    par_list_case_3["log_sigma_log_K_hyper_devs"] <- log(0)
    par_list_case_3["log_sigma_log_L1_hyper_devs"] <- log(0)

    output_data <- simulate_data(obj, par_list_case_3)
  } else if(case == 4) {
    par_list_case_4 <- obj_env_last_par
    par_list_case_4["log_sigma_PE"] <- log(0) # no process error
    par_list_case_4["log_sigma_log_Linf_devs"] <- log(0) # no persistent error
    par_list_case_4["log_sigma_log_K_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_devs"] <- log(0)
    par_list_case_4["log_sigma_log_Linf_hyper_devs"] <- log(0)
    par_list_case_4["log_sigma_log_K_hyper_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_hyper_devs"] <- log(0)

    output_data <- simulate_data(obj, par_list_case_4)

  } else {
    error("provide a value for case that is either 1, 2, 3, or 4")
  }
  return(output_data)
}
