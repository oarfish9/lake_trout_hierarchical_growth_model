# Stebbins 1-15-2024
# helper functions for biphasic_sexspecific_models.R 

# function to write out parameter estimates and standard errors into tibble
format_output_biphasic <- function(sdr, rep) {
  
  fixed <- sdr$par.fixed
  theta_estimates <- unname(sdr$par.fixed)[stringr::str_detect(names(sdr$par.fixed), "theta")] # extract any theta estimates
  
  log_estimates <- unlist(unname(as.list(sdr, what = "Estimate", report = FALSE)[names(fixed)]))
  SEs <- unlist(unname(as.list(sdr, what = "Std. Error", report = FALSE)[names(fixed)]))
  
  estimates <- tibble("parameter" = names(fixed),
                      "log_estimates" = as.numeric(log_estimates),
                      "log_scale_SE" = SEs)
  
  cor_ind <- c(rep$`my_mvn.cov()`[2,1],
               rep$`my_mvn.cov()`[3,1],
               rep$`my_mvn.cov()`[3,2])
  
  cor_ind_names <- c("g_h_ind_cor",
                     "g_L1_ind_cor",
                     "h_L1_ind_cor")
  
  cor_ind_SE <- rep(NA_real_, 3)
  
  correlations_table <- tibble(parameter = cor_ind_names,
                               log_estimates = as.numeric(cor_ind),
                               log_scale_SE = cor_ind_SE)
  
  final_pars <- as_tibble(bind_rows(estimates, correlations_table))
  
  return(final_pars)
}


# create biphasic object for sim
# function for creating obj$sim for MVN_MVN
create_biphasic_obj <- function(sim_vectors) {
  
  # Compile the model and load in the functions
  setwd(here("TMB"))
  if(is.loaded("ltSup_biphasic.cpp")){
    dyn.unload("ltSup_biphasic.cpp")
  }
  TMB::compile("ltSup_biphasic.cpp")
  dyn.load("ltSup_biphasic")
  
  # load the data
  load(here("data-raw", "incdat_2022_biodat.Rdata"))
  idat <- biodat_save
  
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
    log_sigma_log_g_devs = log(0.5),
    log_sigma_log_h_devs = log(0.5), 
    log_sigma_log_L1_devs = log(0.5), 
    log_g_hyper = log(0.2), # 0.2 (Lester 2004), previously 1
    log_h_hyper = log(10),   # 37/42 (Wilson 2019), from Lester 2014 estimate for walleye (0.003), for quince 2008 supplement (50 or 10), h=10 from Lester 2004, mm/year
    log_L1_hyper = log(100),   # starting length for fish at age 1, this is the average of all Ls from idat
    log_g_devs = rep(0, nind), # log of deviations for each individual fish
    log_h_devs = rep(0, nind),
    log_L1_devs = rep(0,nind),
    alpha = 0, # intercept for calculation of age at maturity - previously 83
    log_beta = log(0.1), # slope for calculation of age at maturity - previously 0.5
    theta1 = 0, # set to 0, = no priors; theta1 = g/h
    theta2 = 0, # g/L1
    theta3 = 0, # h/L1
    log_PE = rep(0,nrec), # vector for storing process error, no space for cases where age=0
    log_sigma_PE = log(4) # sigma for process error estimation
  )
  
  # uncomment thetas for a univariate distribution
  map_par <- list(
  )
  
  # Define the random effects
  RE <- c('log_PE', 'log_g_devs', 'log_h_devs', 'log_L1_devs')
  
  # Generate the objective function and gradient
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = "ltSup_biphasic",
    random = RE,
    map = map_par
  )
  
  return(list("obj" = obj))
}

# nested function
simulate_biphasic_data <- function(obj, parameter_list) {
  
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
run_biphasic_cases <- function(obj, obj_env_last_par, case) {
  
  
  if(case == 1) {
    output_data <- simulate_biphasic_data(obj, obj_env_last_par)
    
  } else if(case == 2) {
    par_list_case_2 <- obj_env_last_par
    par_list_case_2["log_sigma_PE"] <- log(0) # no process error
    output_data <- simulate_biphasic_data(obj, par_list_case_2)
    
  } else if(case == 3) {
    par_list_case_3 <- obj_env_last_par
    par_list_case_3["log_sigma_log_g_devs"] <- log(0) # no persistent error
    par_list_case_3["log_sigma_log_h_devs"] <- log(0)
    par_list_case_3["log_sigma_log_L1_devs"] <- log(0)
    
    output_data <- simulate_biphasic_data(obj, par_list_case_3)
  } else if(case == 4) {
    par_list_case_4 <- obj_env_last_par
    par_list_case_4["log_sigma_PE"] <- log(0) # no process error
    par_list_case_4["log_sigma_log_g_devs"] <- log(0) # no persistent error
    par_list_case_4["log_sigma_log_h_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_devs"] <- log(0)
    
    output_data <- simulate_biphasic_data(obj, par_list_case_4)
    
  } else {
    error("provide a value for case that is either 1, 2, 3, or 4")
  }
  return(output_data)
}

# function for running sex-specific model
create_BVN_MVN_obj <- function(sim_vectors, male_dat) {
  
  idat <- male_dat
  
  # Compile the model and load in the functions
  setwd(here("TMB"))
  if(is.loaded("ltSup_MVN_MVN_setnpopsim")){
    dyn.unload("ltSup_MVN_MVN_setnpopsim")
  }
  TMB::compile("ltSup_MVN_MVN_setnpopsim.cpp")
  dyn.load("ltSup_MVN_MVN_setnpopsim")
  
  # Some useful numbers for indexing
  nrec <- nrow(idat)
  nind <- length(unique(idat$idx))
  npop <- length(unique(idat$pop_idx))
  pop_idx_ind <- idat[!duplicated(idat$idx),]$pop_idx #length of nind, values of npop for indexing
  
  # The data to use in the model
  data <- list(
    ind_idx = idat$idx2,
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
  
  map_par <- list(theta5 = factor(NA),
                  theta6 = factor(NA))
  
  # random effects
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
    
    output_data <- simulate_data(obj, par_list_case_3)
  } else if(case == 4) {
    par_list_case_4 <- obj_env_last_par
    par_list_case_4["log_sigma_PE"] <- log(0) # no process error
    par_list_case_4["log_sigma_log_Linf_devs"] <- log(0) # no persistent error
    par_list_case_4["log_sigma_log_K_devs"] <- log(0)
    par_list_case_4["log_sigma_log_L1_devs"] <- log(0)
    
    output_data <- simulate_data(obj, par_list_case_4)
    
  } else {
    error("provide a value for case that is either 1, 2, 3, or 4")
  }
  return(output_data)
}


# old ---------------------------------------------------------------------

# get quantiles of standard error by joining to original model fit
get_relative_error_old <- function(case_x_output_clean, true_estimates_fit) {
  
  cases_not_estimating_thetas <- case_x_output_clean %>% 
    filter(!case %in% c("C", "D", "G", "H"))
  
  theta5_vec <- tibble(parameter = "theta5", log_estimates = as.character(0))
  theta6_vec <- tibble(parameter = "theta6", log_estimates = as.character(0))
  
  estimates_fit_without_thetas <- estimates_fit %>%
    bind_rows(theta5_vec, theta6_vec) %>%
    mutate(true_vals = as.numeric(log_estimates)) %>%
    dplyr::select(true_vals, parameter)
  
  summary_without_thetas <- cases_not_estimating_thetas %>%
    left_join(estimates_fit_without_thetas,
              by = "parameter") %>%
    mutate(log_scale_SE = as.numeric(log_scale_SE),
           log_estimates = as.numeric(log_estimates),
           relative_error = (log_estimates - true_vals) / true_vals,
           median_relative_error = median(relative_error),
           lower_quartile = quantile(relative_error, 0.25),
           upper_quartile = quantile(relative_error, 0.75))
  
  
  if(sum(case_x_output_clean$case %in% c("C", "D", "G", "H")) > 1){
    
    cases_estimating_thetas <- case_x_output_clean %>% 
      filter(case %in% c("C", "D", "G", "H"))
    
    theta5_vec <- tibble(parameter = "theta5", log_estimates = as.character(0.58))
    theta6_vec <- tibble(parameter = "theta6", log_estimates = as.character(0.58))
    
    estimates_fit_with_thetas <- estimates_fit %>%
      bind_rows(theta5_vec, theta6_vec) %>%
      mutate(true_vals = as.numeric(log_estimates)) %>%
      dplyr::select(true_vals, parameter)
    
    summary_with_thetas <- cases_estimating_thetas %>% 
      left_join(estimates_fit_with_thetas,
                by = "parameter") %>%
      mutate(log_scale_SE = as.numeric(log_scale_SE),
             log_estimates = as.numeric(log_estimates),
             relative_error = (log_estimates - true_vals) / true_vals,
             median_relative_error = median(relative_error),
             lower_quartile = quantile(relative_error, 0.25),
             upper_quartile = quantile(relative_error, 0.75))
  } 
  
  
  final_summary <- bind_rows(summary_without_thetas, summary_with_thetas)
  
  return(final_summary)
}

# function to create boxplot
plot_summarized_parallel_output_old <- function(all_cases_with_RE, plot_theta4) {
  if(plot_theta4) {
    p <- all_cases_with_RE %>% 
      filter(!parameter %in% c("K_L1_ind_cor", "Linf_K_ind_cor",
                               "Linf_K_pop_cor", "Linf_L1_ind_cor")) %>% 
      ggplot(aes(x = parameter, y = relative_error)) + geom_boxplot(coef = 0, outlier.shape = NA) + 
      theme_minimal() + 
      facet_wrap(~case, ncol = 1)
  } else {
    p <- all_cases_with_RE %>% 
      filter(!parameter %in% c("theta4", "K_L1_ind_cor", "Linf_K_ind_cor",
                               "Linf_K_pop_cor", "Linf_L1_ind_cor")) %>% 
      ggplot(aes(x = parameter, y = relative_error)) + geom_boxplot(coef = 0, outlier.shape = NA) +
      theme_minimal() + 
      facet_wrap(~case, ncol = 1)
  }
  p + theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    ylab("Relative Error") + ylim(c(-0.25, 0.25))
}

