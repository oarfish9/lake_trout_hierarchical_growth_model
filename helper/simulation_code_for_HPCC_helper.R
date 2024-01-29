# Stebbins 1-15-2024
# helper functions for simulation_code_for_HPCC.R


# function to write out parameter estimates and standard errors into tibble
# argument "MVN_or_BVN" tells the function whether you fit a model estimating
# theta5 and theta6 at the population level (MVN) or not (BVN)
format_output <- function(sdr, rep, MVN_or_BVN) {
  
  fixed <- sdr$par.fixed
  theta_estimates <- unname(fixed)[stringr::str_detect(names(fixed), "theta")] # extract any theta estimates
  
  SEs <- sqrt(diag(sdr$cov.fixed))
  log_estimates <- sdr$par.fixed
  
  
  estimates <- as.data.frame(cbind(names(fixed), log_estimates, SEs))
  colnames(estimates) <- c('parameter', "log_estimates", "log_scale_SE")
  
  # get correlations
  if(MVN_or_BVN == "BVN") {
    cor <- c(rep$`my_mvn.cov()`[2,1],
             rep$`my_mvn.cov()`[3,1],
             rep$`my_mvn.cov()`[3,2],
             rep$`my_mvn_hyper.cov()`[1,2])
    cor_names <- c("Linf_K_ind_cor",
                   "Linf_L1_ind_cor",
                   "K_L1_ind_cor",
                   "Linf_K_pop_cor")
    cor_SEs <- rep(NA, 4)
  } else if(MVN_or_BVN == "MVN") {
    cor <- c(rep$`my_mvn.cov()`[2,1],
             rep$`my_mvn.cov()`[3,1],
             rep$`my_mvn.cov()`[3,2],
             rep$`my_mvn_hyper.cov()`[1,2],
             rep$`my_mvn_hyper.cov()`[1,3],
             rep$`my_mvn_hyper.cov()`[2,3])
    cor_names <- c("Linf_K_ind_cor",
                   "Linf_L1_ind_cor",
                   "K_L1_ind_cor",
                   "Linf_K_pop_cor",
                   "Linf_L1_pop_cor",
                   "L1_K_pop_cor")
    cor_SEs <- rep(NA, 6)
  }
  
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

# function for creating obj$sim for MVN_MVN
# this creates an object "obj" that expects a theta5 and theta6 value 
# and gets called to simulate datasets within simulations
create_MVNMVN_obj <- function(idat, sim_vectors) {
  
  # Compile the model and load in the functions
  setwd(here::here("TMB"))
  if(is.loaded("ltSup_MVN_MVN_setnpopsim.cpp")){
    dyn.unload("ltSup_MVN_MVN_setnpopsim.cpp")
  }
  TMB::compile("ltSup_MVN_MVN_setnpopsim.cpp")
  dyn.load("ltSup_MVN_MVN_setnpopsim")
  
  # load the data
  load(here::here("data", "laketrout_inc_dat.Rdata"))
  
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
    log_sigma = log(0.004),
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
    log_sigma_PE = log(0.2), # sigma for process error estimation
    log_Linf_devs = rep(0, nind), # log of deviations for each individual fish
    log_K_devs = rep(0, nind),
    log_L1_devs = rep(0, nind),
    log_Linf_hyper_devs = rep(0, npop), # log of deviations for population dist.
    log_K_hyper_devs = rep(0, npop),
    log_L1_hyper_devs = rep(0, npop),
    log_PE = rep(0, nrec) # vector for storing process error, no space for cases where age=0
  )
  
  # uncomment thetas for a univariate distribution
  map_par <- list(
  )
  
  # Define the random effects
  RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
          'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')
  
  # Generate the objective function and gradient
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = 'ltSup_MVN_MVN_setnpopsim',
    random = RE,
    map = map_par
  )
  
  return(list("obj" = obj))
}

# this takes in args for npop_sim, fixed, theta5, and theta6 to run 
# a simulation (simulate data and fit) according to specifications
simulate_and_fit <- function(original_parameter_estimates, sim_vectors,
                             theta5, theta6, npop_sim, maxage,
                             fixed, obj_sim){
  
  # change the parameters used to simulate data to have the inputted theta5 or theta6 values
  pars_sim <- original_parameter_estimates
  pars_sim[["theta5"]] <- theta5
  pars_sim[["theta6"]] <- theta6
  
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
  
  # use good parameter estimates as starting values
  pars_sim_new <- list(
    log_sigma = unname(pars_sim["log_sigma"]),
    log_sigma_log_Linf_devs = unname(pars_sim["log_sigma_log_Linf_devs"]),
    log_sigma_log_K_devs = unname(pars_sim["log_sigma_log_K_devs"]), 
    log_sigma_log_L1_devs = unname(pars_sim["log_sigma_log_L1_devs"]),
    log_sigma_log_Linf_hyper_devs = unname(pars_sim["log_sigma_log_Linf_hyper_devs"]),
    log_sigma_log_K_hyper_devs = unname(pars_sim["log_sigma_log_K_hyper_devs"]),
    log_sigma_log_L1_hyper_devs = unname(pars_sim["log_sigma_log_L1_hyper_devs"]),
    log_Linf_hyper = unname(pars_sim["log_Linf_hyper"]),
    log_K_hyper = unname(pars_sim["log_K_hyper"]),
    log_L1_hyper = unname(pars_sim["log_L1_hyper"]),
    theta1 = unname(pars_sim["theta1"]),
    theta2 = unname(pars_sim["theta2"]),
    theta3 = unname(pars_sim["theta3"]),
    theta4 = unname(pars_sim["theta4"]),
    theta5 = unname(pars_sim["theta5"]),
    theta6 = unname(pars_sim["theta6"]),
    log_sigma_PE = unname(pars_sim["log_sigma_PE"]),
    log_Linf_devs = rep(0, sim_vectors$nind_sim),
    log_K_devs = rep(0, sim_vectors$nind_sim),
    log_L1_devs = rep(0, sim_vectors$nind_sim),
    log_Linf_hyper_devs = rep(0, npop_sim),
    log_K_hyper_devs = rep(0, npop_sim),
    log_L1_hyper_devs = rep(0, npop_sim),
    log_PE = rep(0, sim_vectors$nrec_sim)
  )
  
  if(fixed == TRUE){
    # don't estimate theta5 and theta6
    map_par <- list(
      theta5 = factor(NA),
      theta6 = factor(NA)
    )
  } else {
    map_par <- list()
  }
  
  RE <- c('log_Linf_devs', 'log_K_devs','log_L1_devs',
          'log_Linf_hyper_devs', 'log_K_hyper_devs', 'log_L1_hyper_devs', 'log_PE')
  
  obj_sim <- MakeADFun(
    data = data_sim,
    parameters = pars_sim_new,
    DLL = 'ltSup_MVN_MVN_setnpopsim',
    random = RE,
    map = map_par
  )
  
  # if model fit errors out, set error_out to TRUE
  error_out <- tryCatch({nlminb(start = obj_sim$par, 
                                objective = obj_sim$fn, 
                                gradient = obj_sim$gr)
    error_out <- FALSE},
    error = function(e) {
      error_out <- TRUE
    })
  
  # if model errored out, return an tibble filled with NAs
  if(error_out) {
    estimates_sim <- tibble("parameter" = rep(NA, 20),
                            "log_estimates" = rep(NA, 20),
                            "log_scale_SE" = rep(NA, 20))
  } else{
    rep_sim <- obj_sim$report()
    sdr_sim <- sdreport(obj_sim)
    
    # output
    if(fixed == TRUE){
      estimates_sim <- format_output(sdr_sim, rep_sim, MVN_or_BVN = "BVN")
    } else {
      estimates_sim <- format_output(sdr_sim, rep_sim, MVN_or_BVN = "MVN")
    }
  }
  
  
  return(estimates_sim)
}


# this calls simulate_and_fit and run_simulations according to the 
# specified number of simulations (nsims)
run_simulations <- function(nsims, original_parameter_estimates,
                            sim_vectors, theta5, theta6,
                            npop_sim, maxage, fixed, obj_sim){
  output <- list()
  converged <- as.logical(rep(NA, nsims))
  
  for(i in 1:nsims) {
    
    output[[i]] <- simulate_and_fit(original_parameter_estimates,
                                    sim_vectors, theta5,
                                    theta6, npop_sim, maxage,
                                    fixed = fixed,
                                    obj_sim) |> 
      dplyr::mutate(sim = i,
                    log_scale_SE = ifelse(log_scale_SE == "NaN", NA, log_scale_SE))
    
    # check if simulate_and_fit returned an empty tibble (error out)
    if(sum(is.na(output[[i]]$parameter)) > 0) {
      converged[i] <- FALSE
      # check if simulate_and_fit returned a table with more than one NA for
      # log_scale_SE
    } else if(fixed == TRUE & sum(is.na(output[[i]]$log_scale_SE)) > 4) {
      converged[i] <- FALSE
    } else if(fixed == FALSE & sum(is.na(output[[i]]$log_scale_SE)) > 6) {
      converged[i] <- FALSE
    } else{
      converged[i] <- TRUE
    }
  }
  
  print(converged)
  raw_output <- do.call(rbind, output)
  
  return(list("raw_output" = raw_output,
              "prop_converged" = sum(converged)/nsims))
  
}