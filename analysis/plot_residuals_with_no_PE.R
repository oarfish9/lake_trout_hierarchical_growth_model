# Stebbins 2-23-2024
# This takes the best fitting model and produces a predicted length with no
# process error effect (for residual plotting)

# libraries
library(tidyverse)
library(here)
library(TMB)

# helper functions
source(here("helper", "best_fitting_model_helper.R"))

# Compile the model and load in the functions
setwd(here("TMB"))
if(is.loaded("ltSup_MVN_BVN_v2_no_PE")){
  dyn.unload("ltSup_MVN_BVN_v2_no_PE")
}
TMB::compile("ltSup_MVN_BVN_v2_no_PE.cpp")
dyn.load("ltSup_MVN_BVN_v2_no_PE")

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
  DLL = 'ltSup_MVN_BVN_v2_no_PE',
  random = RE,
  map = map_par
)

# run the model
opt <- nlminb(start=obj$par, objective=obj$fn, 
              gradient=obj$gr)

# report out
rep <- obj$report()
sdr <- sdreport(obj)

# get residuals
load(here("data", "fish_ids_with_centroids.Rdata"))

residuals <- tibble("id" = idat$id,
                   "pop_id" = idat$pop_idx + 1,
                   "residuals" = rep$resid,
                   "log_PE" = log(rep$PE),
                   "pred_L_no_PE" = rep$L_hat_no_PE,
                   "delta" = rep$delta,
                   "pred_L_w_PE" = rep$L_hat,
                   "log_residuals_no_PE" = rep$log_resid_no_PE,
                   #"log_residuals_no_PE" = log(rep$L) - log(rep$L_hat_no_PE),
                   "log_residuals" = log(rep$L) - log(rep$L_hat),
                   "age" = idat$Age) |> 
  left_join(fish_ids_with_centroids |> 
              select(id, location))

# plot
top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")

log_resid_plot_no_PE <- residuals |> 
  rename(Location = location) |> 
  ggplot(aes(x = age, log_residuals_no_PE, color = Location)) + 
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Location) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = top_palette) +
  theme_bw() +
  labs(x = "Age", y = "Log-scale residuals",
       title = "Log-scale residuals of predicted (less process error) vs. observed (back-calculated) length by population")

ggsave(log_resid_plot_no_PE, file = here("figures", "log_resid_plot_no_PE.png"))

PE_by_age_plot <- residuals |> 
  rename(Location = location) |> 
  ggplot(aes(x = age, log_PE, color = Location)) + 
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~Location) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Estimated log-scale process errors")

ggsave(PE_by_age_plot, file = here("figures", "PE_by_age_plot.png"))

diff_between_methods <- residuals |> 
  mutate(diff_between_methods = pred_L_w_PE - pred_L_no_PE) |> 
  rename(Location = location) |> 
  ggplot(aes(x = age, diff_between_methods, color = Location)) + 
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~Location) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Difference (mm)")

ggsave(diff_between_methods, file = here("figures", "diff_between_methods.png"))
  

pred_L_by_PE_plot <- residuals |> 
  rename(Location = location) |> 
  pivot_longer(c(pred_L_no_PE, pred_L_w_PE),
                 names_to = "Method",
                 values_to = "pred_L") |> 
  mutate(Method = ifelse(Method == "pred_L_no_PE", 
                         "Predicted length (no process error)",
                         "Predicted length (with process error)")) |> 
  mutate(group_id = paste(id, Method, sep = "_")) |> 
  ggplot(aes(x = age, pred_L, group = group_id, color = Method)) + 
  geom_line(alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = top_palette[2:3]) +
  facet_wrap(~Location, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age",
       y = "Predicted length (mm)")

ggsave(pred_L_by_PE_plot, file = here("figures", "pred_L_by_PE_plot.png"))

pred_delta <- residuals |> 
  rename(Location = location) |> 
  # ggplot(aes(x = delta, fill = Location)) +
  # geom_histogram() +
  ggplot(aes(x = age, y = delta, color = Location)) +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~Location, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age",
       y = "Predicted length increment (delta)") +
  scale_color_manual(values = top_palette)

ggsave(pred_delta, file = here("figures", "pred_delta.png"))

