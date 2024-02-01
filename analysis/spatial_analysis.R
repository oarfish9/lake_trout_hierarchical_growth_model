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
                   "residuals" = rep$resid) |> 
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

plot_inds_with_distance |> 
  pivot_longer(Linf:L1, names_to = "param", values_to = "estimate") |> 
  ggplot(aes(x = location, y = estimate, fill = location)) +
  geom_boxplot() +
  facet_wrap(~param, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_manual(values = top_palette)

# plot residuals and color by population
# are these fitted vB trajectories vs back-calculated?

residuals |> 
  ggplot(aes(x = id, residuals, color = location)) + 
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~location, scales = "free") +
  xlab("Fish ID") +
  ggtitle("Residuals of predicted vs. observed (back-calculated) length by population") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
