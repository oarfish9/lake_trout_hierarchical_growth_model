# Stebbins 1-15-2024
# this contains helper functions for running the script 
# evaluate_alternate_models.R

# functions ---------------------------------------------------------------

# function to write out parameter estimates and standard errors into tibble
format_output <- function(sdr, rep, 
                          ind_level = c("UVN", "BVN", "MVN"), 
                          pop_level = c("UVN", "BVN", "MVN")) {
  
  fixed <- sdr$par.fixed # get fixed pars from sdr (estimates)
  theta_estimates <- unname(sdr$par.fixed)[stringr::str_detect(names(sdr$par.fixed), "theta")] # extract any theta estimates
  # extract fixed parameter SEs from sdrsim report
  # https://www.rdocumentation.org/packages/TMB/versions/1.9.0/topics/as.list.sdreport
  # SEs <- unlist(unname(as.list(sdr, what = "Std. Error", report = FALSE)[1:length(fixed)]))
  # log_estimates <- unlist(unname(as.list(sdr, what = "Estimate", report = FALSE)[1:length(fixed)]))
  
  log_estimates <- unlist(unname(as.list(sdr, what = "Estimate", report = FALSE)[names(fixed)]))
  SEs <- unlist(unname(as.list(sdr, what = "Std. Error", report = FALSE)[names(fixed)]))
  
  estimates <- as.data.frame(cbind(names(fixed), log_estimates, SEs))
  #estimates <- as.data.frame(cbind(names(fixed), unname(fixed), SEs))
  colnames(estimates) <- c('parameter', "log_estimates", "log_scale_SE")
  
  # get correlations
  if(ind_level == "UVN") {
    cor_ind <- NULL
    cor_ind_names <- NULL
    cor_ind_SE <- NULL
  } else if(ind_level == "BVN") {
    cor_ind <- c(rep$`my_mvn.cov()`[2,1])
    cor_ind_names <- c("Linf_K_ind_cor")
    cor_ind_SE <- rep(NA, 1)
  } else if(ind_level == "MVN") {
    cor_ind <- c(rep$`my_mvn.cov()`[2,1],
                 rep$`my_mvn.cov()`[3,1],
                 rep$`my_mvn.cov()`[3,2])
    cor_ind_names <- c("Linf_K_ind_cor",
                       "Linf_L1_ind_cor",
                       "K_L1_ind_cor")
    cor_ind_SE <- rep(NA, 3)
  }
  
  
  if(pop_level == "UVN") {
    cor_pop <- NULL
    cor_pop_names <- NULL
    cor_pop_SE <- NULL
  } else if(pop_level == "BVN") {
    cor_pop <- c(rep$`my_mvn_hyper.cov()`[2,1])
    cor_pop_names <- c("Linf_K_pop_cor")
    cor_pop_SE <- rep(NA, 1)
  } else if(pop_level == "MVN") {
    cor_pop <- c(rep$`my_mvn_hyper.cov()`[2,1],
                 rep$`my_mvn_hyper.cov()`[3,1],
                 rep$`my_mvn_hyper.cov()`[3,2])
    cor_pop_names <- c("Linf_K_pop_cor",
                       "Linf_L1_pop_cor",
                       "K_L1_pop_cor")
    cor_pop_SE <- rep(NA, 3)
  }
  
  cor <- c(cor_ind, cor_pop)
  cor_names <- c(cor_ind_names, cor_pop_names)
  cor_SEs <- c(cor_ind_SE, cor_pop_SE)
  
  
  correlations_table <- as_tibble(cbind(cor_names, cor, cor_SEs))
  names(correlations_table) <- c("parameter", "log_estimates", "log_scale_SE")
  
  final_pars <- as_tibble(bind_rows(estimates, correlations_table))
  
  return(final_pars)
}

# function to convert thetas to correlation values
thetas_to_cor <- function(thetas) {
  # input is a vector of thetas (right now this only works when you have 3 thetas)
  # equation found here: https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html
  
  # create L matrix
  L <- matrix(0, nrow=3, ncol=3)
  diag(L) <- rep(1, 3)
  lower_tri <- lower.tri(L)
  t_lower_tri <- t(lower_tri)
  L[t_lower_tri] <- thetas
  
  # create D
  middle <- L %*% t(L)
  D <- diag(middle)
  Dinvsqrt <- diag(1/sqrt(D))
  
  # create correlation matrix
  cormat <- Dinvsqrt %*% middle %*% Dinvsqrt
  return(cormat)
}

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



# this takes in args for npop_sim, fixed, theta5, and theta6 to run 
# a simulation (simulate data and fit) according to specifications
simulate_and_fit <- function(original_parameter_estimates, sim_vectors,
                             theta5, theta6, npop_sim, maxage,
                             fixed, obj_sim){
  
  # change the parameters used to simulate data to have the inputted theta5 or theta6 values
  pars_sim <- original_parameter_estimates
  pars_sim[[15]] <- theta5
  pars_sim[[16]] <- theta6
  
  # simulate a new dataset
  sim_dat <- obj_sim$simulate(pars_sim, complete = T)
  
  # fit the new dataset:
  
  # set new data and parameters
  data_sim <- list(
    ind_idx = sim_vectors$ind_idx_sim, 
    pop_idx_ind = sim_vectors$pop_idx_ind_sim, 
    ind_idx_sim = sim_vectors$ind_idx_sim,
    pop_idx_ind_sim = sim_vectors$pop_idx_ind_sim,
    npop = npop_sim,
    npop_sim = npop_sim,
    nrec = sim_vectors$nrec_sim, 
    nind = maxage * npop_sim,
    L = sim_dat$sim_L_hat, # simulated L
    age = sim_dat$age_sim,
    maxage = maxage,
    age_sim = sim_vectors$age
  )
  
  pars_sim <- list(
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
    log_Linf_devs = rep(0, sim_vectors$nind_sim), # log of deviations for each individual fish
    log_K_devs = rep(0, sim_vectors$nind_sim),
    log_L1_devs = rep(0, sim_vectors$nind_sim),
    log_Linf_hyper_devs = rep(0, npop_sim), # log of deviations for population dist.
    log_K_hyper_devs = rep(0, npop_sim),
    log_L1_hyper_devs = rep(0, npop_sim),
    log_PE = rep(0, sim_vectors$nrec_sim) # vector for storing process error, no space for cases where age=0
  )
  
  if(fixed == TRUE){
    # don't estimate theta5 and theta6
    map_par <- list(
      theta5 = factor(NA),
      theta6 = factor(NA)
    )
  } else {
    map_par <- list() # estimate theta5 and theta6
  }
  
  obj_sim <- MakeADFun(
    data = data_sim,
    parameters = pars_sim,
    DLL = 'ltSup_MVN_MVN_setnpopsim',
    random = RE,
    map = map_par
  )
  
  # if model fit errors out, keep going
  # of a tibble
  tryCatch(optsimtest <- nlminb(start = obj_sim$par, 
                                objective = obj_sim$fn, 
                                gradient = obj_sim$gr),
           error = function(e) {return("model did not converge")})
  
  
  rep_sim <- obj_sim$report()
  sdr_sim <- sdreport(obj_sim)
  
  # output
  if(fixed == TRUE){
    estimates_sim <- format_output(sdr_sim, rep_sim, MVN_or_BVN = "BVN")
  } else {
    estimates_sim <- format_output(sdr_sim, rep_sim, MVN_or_BVN = "MVN")
  }
  return(estimates_sim)
}


# get quantiles of standard error by joining to original model fit
get_relative_error <- function(case_x_output_clean, true_estimates_fit) {
  
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

# function to format tibble with exp of estimates -------------------------

get_estimates <- function(model_results) {
  results <- model_results$estimates_fit |> 
    mutate(estimates = case_when(str_detect(log_estimates, "cor") ~ log_estimates,
                                 TRUE ~ exp(log_estimates)),
           AIC = model_results$AIC,
           convergence = case_when(!is.na(log_estimates) & is.na(log_scale_SE) ~ "NO",
                                   TRUE ~ "YES"))
  return(results)
  
}
