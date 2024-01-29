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