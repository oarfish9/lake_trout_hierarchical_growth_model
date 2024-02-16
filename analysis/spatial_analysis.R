# spatial analysis

library(tidyverse)
library(ggridges)

load(here("data", "best_fitting_rep.Rdata"))
load(here("data-raw", "incdat_2021.Rdata"))
load(here("data", "fish_ids_with_centroids.Rdata"))
load(here("data", "lake_superior_lat_longs.Rdata"))
load(here("data", "nets_with_lat_longs.Rdata"))
pop_ids <- read_csv(here("data", "Lake Trout Lake Names and IDs.csv"))

top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")


ind_Linfs <- c()
ind_Ks <- c()
ind_L1s <- c()
individual_pops <- idat |> distinct(id, pop_idx) |> pull(pop_idx)
individual_pops <- individual_pops + 1
for(i in 1:length(rep$log_Linf_devs)) {
  ind_Linfs[i] <- exp(rep$log_Linf_hyper + rep$log_Linf_hyper_devs[individual_pops[i]] + rep$log_Linf_devs[i])
  ind_Ks[i] <- exp(rep$log_K_hyper + rep$log_K_hyper_devs[individual_pops[i]] + rep$log_K_devs[i])
  ind_L1s[i] <- exp(rep$log_L1_hyper + rep$log_L1_hyper_devs[individual_pops[i]] + rep$log_L1_devs[i])
}

ind_fish <- idat |> 
  distinct(id) |> 
  left_join(fish_ids_with_centroids, by = "id")

# residuals
residuals = tibble("id" = idat$id,
                   "pop_id" = idat$pop_idx + 1,
                   "residuals" = rep$resid,
                   "log_residuals" = log(rep$L) - log(rep$L_hat),
                   "age" = idat$Age) |> 
  left_join(fish_ids_with_centroids |> 
              select(id, location))


# plots -------------------------------------------------------------------
plot_inds_with_distance <- tibble("Linf" = ind_Linfs,
                                  "K" = ind_Ks,
                                  "L1" = ind_L1s,
                                  "id" = ind_fish$id,
                                  "lat" = ind_fish$latitude,
                                  "long" = ind_fish$longitude,
                                  "location" = ind_fish$location,
                                  "year_capture" = ind_fish$year_capture)

plot_inds_with_distance |> 
  ggplot(aes(x = Linf, y = K, color = location)) + 
  geom_point()

plot_inds_with_distance |> 
  ggplot(aes(x = Linf , y = location,  fill = location)) +
  geom_density_ridges(scale = 1)

vb_estimates_by_loc <- plot_inds_with_distance |> 
  rename(Location = location) |> 
  pivot_longer(Linf:L1, names_to = "param", values_to = "estimate") |> 
  mutate(symbols = case_when(param == "Linf" ~ "L[infinity]",
                             param == "L1" ~ "L[1]",
                            TRUE ~ param)) |> 
  ggplot(aes(x = Location, y = estimate, fill = Location)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~symbols, scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = top_palette) +
  labs(y = "Estimated parameter value",
       title = "Distribution of von Bertalanffy parameter estimates by location")

ggsave(vb_estimates_by_loc, file = here("figures", "vb_estimates_by_loc.png"))

# plot residuals and color by population
# plot the x axis as "distance between locations" ?

log_resid_plot <- residuals |> 
  rename(Location = location) |> 
  ggplot(aes(x = age, log_residuals, color = Location)) + 
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-0.001, 0.001)) +
  facet_wrap(~Location) +
  theme(legend.position = "below") +
  scale_color_manual(values = top_palette) +
  theme_bw() +
  labs(x = "Age", y = "Log-scale residuals (mm)",
       title = "Log-scale residuals of predicted vs. observed (back-calculated) length by population")

ggsave(log_resid_plot, file = here("figures", "log_resid_plot.png"))

residuals |> 
  filter(residuals > -0.01 & residuals <0.01) |> 
  ggplot(aes(x = residuals , y = location,  fill = location)) +
  geom_density_ridges(scale = 1) +
  ggtitle("Residuals by location")
