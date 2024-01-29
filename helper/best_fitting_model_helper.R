# helper functions for best_fitting_model.R

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