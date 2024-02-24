# create plots for figures

# libraries
library(viridis)
library(tidyverse)
library(latex2exp)
library(here)

# set palette
top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")
#top_palette <- inferno(20)

backup_palette <- turbo(10)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load files --------------------------------------------------------------
load(here("data", "best_fitting_model_with_PIs.RData"))
load(here("data", "HPCC_relative_error_with_CIs_v2.Rdata"))
load(here("data", "prop_converged_v2.Rdata"))
load(here("data", "persistent_transient_sim_results.Rdata"))
load(here("data", "persistent_transient_sim_results_biphasic.Rdata"))
load(here("data", "persistent_transient_sim_results_sex_specific.Rdata"))


# plot best fitting model params ------------------------------------------
confidence_intervals_plot <- probability_intervals |> 
  filter(parameter %in% c("log_sigma_log_Linf_devs",
                          "log_sigma_log_L1_devs",
                          "log_sigma_log_K_devs",
                          "log_sigma_log_Linf_hyper_devs",
                          "log_sigma_log_L1_hyper_devs",
                          "log_sigma_log_K_hyper_devs",
                          "log_Linf_hyper",
                          "log_L1_hyper",
                          "log_K_hyper")) |> 
  select(parameter, log_estimates, estimates) 

CI_plot_vals <- tibble(parameter = rep(c("Linf", "K", "L1"), 1, each = 2),
                       level = rep(c("Individual", "Population"), 3),
                       log_est = rep(c(confidence_intervals_plot$log_estimates[7],
                                       confidence_intervals_plot$log_estimates[8],
                                       confidence_intervals_plot$log_estimates[9]), 1, each = 2),
                       sds = c(confidence_intervals_plot$log_estimates[1],
                               confidence_intervals_plot$log_estimates[4],
                               confidence_intervals_plot$log_estimates[2],
                               confidence_intervals_plot$log_estimates[5],
                               confidence_intervals_plot$log_estimates[3],
                               confidence_intervals_plot$log_estimates[6]),
                       symbols = rep(c("L[infinity]", "K", "L[1]"), 1, each = 2)) |> 
  mutate(sds = exp(sds),
         est = exp(log_est),
         upper_CI = exp(log_est + (2 * sds)),
         lower_CI = exp(log_est - (2 * sds))) |> 
  glimpse()

ind_pop_plot <- ggplot(CI_plot_vals, aes(x = level, y = est, color = level)) +
  geom_errorbar(aes(ymin = lower_CI, 
                    ymax = upper_CI),
                width = 0.1, 
                linewidth = 0.8) + 
  geom_point() +
  scale_color_manual(values = top_palette[2:3]) +
  facet_wrap(~symbols, scales = "free", labeller = label_parsed) +
  theme_bw() +
  ylab("Median and 95% probability intervals") +
  xlab("Parameter and level of variation")

ggsave(ind_pop_plot, file = here("figures", "ind_pop_plot.png"))

# plot HPCC results -------------------------------------------------------

# plotting function for thetas
plot_relative_errors <- function(relative_error_CIs, plot_thetas) {
  if(plot_thetas) {
    new_dat <- relative_error_CIs |> 
      filter(!str_detect(parameter, "cor")) |> 
      filter(!parameter %in% c("theta1", "theta2", "theta3",
                               "theta4", "theta5", "theta6"))
    
    p <- ggplot(new_dat, aes(x = parameter, median)) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), 
                    width = 0.1, 
                    size = 0.8) +
      geom_point() + 
      facet_wrap(~case, nrow = 4) + 
      theme(strip.background = element_blank()) +
      theme_bw() +
      scale_x_discrete(labels = c(
        expression(K[...]),
        expression(L1[...]),
        expression(L[infinity[...]]),
        expression(sigma[epsilon]),
        expression(sigma[K[i]]),
        expression(sigma[K[p]]),
        expression(sigma[L1[i]]),
        expression(sigma[L1[p]]),
        expression(sigma[L[infinity[i]]]),
        expression(sigma[L[infinity[p]]]),
        expression(sigma[delta]),
        expression(theta[1]),
        expression(theta[2]),
        expression(theta[3]),
        expression(theta[4]),
        expression(theta[5]),
        expression(theta[6]))) +
      ggtitle("Distribution of relative error by case, non-correlation parameters")
  } else {
    new_dat <- relative_error_CIs |>  
      filter(!str_detect(parameter, "cor")) |> 
      filter(parameter %in% c("theta1", "theta2", "theta3",
                              "theta4", "theta5", "theta6"))
    
    p <- ggplot(new_dat, aes(x = parameter, median)) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), 
                    width = 0.1, 
                    size = 0.8) +
      geom_point() + 
      facet_wrap(~case, scales = "free", nrow = 4) + 
      theme(strip.background = element_blank()) +
      theme_bw() +
      scale_x_discrete(labels = c(
        expression(theta[1]),
        expression(theta[2]),
        expression(theta[3]),
        expression(theta[4]),
        expression(theta[5]),
        expression(theta[6]))) +
      ggtitle("Distribution of relative error by case, correlation parameters")
  }
  p + 
    ylab("Relative Error") +
    xlab("Parameter") +
    theme(text=element_text(size=12))
}

theta_plot <- plot_relative_errors(relative_error_CIs, FALSE)
non_theta_plot <- plot_relative_errors(relative_error_CIs, TRUE)

ggsave(theta_plot, file = here("figures", "relative_errors_thetas.png"))
ggsave(non_theta_plot, file = here("figures", "relative_errors_non_thetas.png"))

# plotting function for correlations
# plotting function
plot_relative_errors_cor <- function(relative_error_CIs, plot_cors) {
  if(plot_cors == FALSE) {
    new_dat <- relative_error_CIs |> 
      filter(!str_detect(parameter, "theta")) |> 
      filter(!parameter %in% c("K_L1_ind_cor", "L1_K_pop_cor", "Linf_K_ind_cor",
                               "Linf_K_pop_cor", "Linf_L1_ind_cor", "Linf_L1_pop_cor"))
    
    p <- ggplot(new_dat, aes(x = parameter, median)) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), 
                    width = 0.1, 
                    size = 0.8) +
      geom_point() + 
      facet_wrap(~case, nrow = 4) + 
      theme(strip.background = element_blank()) +
      theme_bw() +
      scale_x_discrete(labels = c(
        expression(K[...]),
        expression(L[1[...]]),
        expression(L[infinity[...]]),
        expression(sigma[epsilon]),
        expression(sigma[K[i]]),
        expression(sigma[K[p]]),
        expression(sigma[L[1[i]]]),
        expression(sigma[L[1[p]]]),
        expression(sigma[L[infinity[i]]]),
        expression(sigma[L[infinity[p]]]),
        expression(sigma[delta])))
  } else if(plot_cors == TRUE) {
    
    new_dat <- relative_error_CIs |> 
      filter(parameter %in% c("K_L1_ind_cor", "L1_K_pop_cor", "Linf_K_ind_cor",
                              "Linf_K_pop_cor", "Linf_L1_ind_cor", "Linf_L1_pop_cor"))
    
    p <- ggplot(new_dat, aes(x = parameter, median)) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), 
                    width = 0.1, 
                    size = 0.8) +
      geom_point() + 
      facet_wrap(~case, scales = "free_y", nrow = 4) +
      scale_x_discrete(labels = c(
        expression(L[infinity[i]]/K[i]),
        expression(L[infinity[i]]/L[1[i]]),
        expression(L[1[i]]/K[i]),
        expression(L[infinity[p]]/K[p]),
        expression(L[infinity[p]]/L[1[p]]),
        expression(L[1[p]]/K[p]))) +
      theme(strip.background = element_blank()) +
      theme_bw()
  }
  p + 
    ylab("Relative Error") +
    xlab("Parameter") +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

non_cor_plot <- plot_relative_errors_cor(relative_error_CIs, FALSE)
cor_plot <- plot_relative_errors_cor(relative_error_CIs, TRUE)

ggsave(cor_plot, file = here("figures", "relative_errors_cors.png"),
       width = 7, height = 6)
ggsave(non_cor_plot, file = here("figures", "relative_errors_non_cors.png"),
       width = 7, height = 5)


# plot pers vs trans ------------------------------------------------------

pers_trans <- function(all_cases) {
  all_cases |> 
    filter(age %in% c(4, 15, 40),
           sim_length < 1000,
           case != "4") |> 
    mutate(case = case_when(case == "1" ~ "Per + Trans",
                            case == "2" ~ "Per",
                            case == "3" ~ "Trans",
                            TRUE ~ case),
           age_title = paste0("Age = ", age),
           age_title_f = factor(age_title, levels = c("Age = 4", "Age = 15", "Age = 40"))) |> 
    ggplot(aes(x = sim_length, color = case)) + 
    geom_density(size = 1.2) +
    facet_wrap(~age_title_f, nrow = 3, scales = "free") +
    scale_color_manual(values = top_palette) + 
    #theme_minimal() +
    theme_bw() +
    #ggtitle(paste0("Density plot of true lengths for age by variation simulation case")) +
    xlab('Simulated True Length') +
    ylab("Density") +
    labs(color = "Case")
}

facet_plot <- pers_trans(all_cases)
ggsave(facet_plot, file = here("figures", "pers_trans_plot.png"))

biphasic_facet_plot <- pers_trans(all_cases_biphasic)
ggsave(biphasic_facet_plot, file = here("figures", "pers_trans_plot_biphasic.png"))

sex_specific_facet_plot <- pers_trans(all_cases_sex_specific)
ggsave(sex_specific_facet_plot, file = here("figures", "pers_trans_plot_sex_specific.png"))

# plot prop converged -----------------------------------------------------

prop_converged_plot <- prop_converged |> 
  group_by(case) |> 
  summarise(median = median(prop_converged, na.rm = T),
            lower = quantile(prop_converged, 0.05),
            upper = quantile(prop_converged, 0.95)) |> 
  ggplot(aes(x = case, y = median)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.1, 
                size = 0.8) +
  theme_bw() +
  geom_point() +
  xlab("Simulation case") +
  ylab("Proportion converged") +
  ggtitle("Proportion of simulations converged (median and 5%/95% CIs)")
ggsave(prop_converged_plot, file = here("figures", "prop_converged_plot.png"))


# scratch -----------------------------------------------------------------

ggplot(df_params_L1, aes(x=params_L1, y=exp_means_L1)) + 
  theme_bw() + labs(x='', y='L1 (mm)') + 
  geom_errorbar(aes(ymin=exp(log_means_L1-(1.96*sds_L1)), 
                    ymax=exp(log_means_L1+(1.96*sds_L1))),
                color=c('red', 'orange'), width=0.1, size=0.8) + 
  geom_point()



# plot function
plot_pers_vs_trans <- function(data, age_arg) {
  data |> 
    filter(age == age_arg) |> 
    filter(sim_length < 1000) |> 
    ggplot(aes(x = sim_length, color = case)) + 
    geom_density(size = 1.2) + 
    facet_wrap(~age_arg) +
    scale_color_manual(values = top_palette) + 
    theme_minimal() +
    ggtitle(paste0("Density plot of true lengths for age ", age_arg, 
                   " by variation simulation case")) +
    xlab('Simulated True Length')
}

age_4_plot <- all_cases |> 
  filter(case != "4") |>
  plot_pers_vs_trans(age = 4)

age_15_plot <- all_cases |> 
  filter(case != 4) |>
  plot_pers_vs_trans(age = 15)

age_40_plot <- all_cases |> 
  filter(case != 4) |>
  plot_pers_vs_trans(age = 40)

ggsave(age_4_plot, file = here("figures", "age_4_plot.png"))
ggsave(age_15_plot, file = here("figures", "age_15_plot.png"))
ggsave(age_40_plot, file = here("figures", "age_40_plot.png"))


# map plots ---------------------------------------------------------------

load("data/lake_superior_lat_longs.Rdata")
load("data/nets_with_lat_longs.Rdata")

lake_superior_plot <- lake_superior_lat_longs |>
  ggplot(aes(x = longitude, y = latitude, group = location)) +
  geom_polygon(fill = 'white', colour = "grey50") + 
  geom_polygon(data = lake_superior_lat_longs |> 
                 filter(location == "Isle_Royale"), 
               aes(x = longitude, y = latitude, group = location),
               fill = "white", color = "grey50") +
  geom_point(data = nets_with_lat_longs, 
             aes(x = longitude, y = latitude, color = location),
             alpha = 0.6, size = 2, stroke = 1) +
  scale_color_manual(values = top_palette, name = "Population") +
  labs(x = "", y = "", 
       title = "Lake Superior Sampling Locations") +
  theme(legend.position = "bottom")

ggsave(lake_superior_plot, file = here("figures", "lake_superior_plot.png"))

