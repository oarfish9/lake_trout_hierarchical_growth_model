# back-calculation methods

rm(list=ls())
# libraries
library(dplyr)
library(ggplot2)
library(here)
library(wesanderson)
library(smatr) # standardized major axis regression


# plot otolith v body -----------------------------------------------------
load(here("data-raw", "incdat_2021.Rdata"))
short_palette = c("#354823", "#D69C4E")
top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")

backup_palette <- turbo(10)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# we don't have the individual trajectories; we can only show
# fit of our linear model to the data

total_length_radius_linear <- lm(TL_mm ~ M_radius, data = idat)
standardized_major_axis <- with(idat, line.cis(y = TL_mm, x = M_radius, method = "SMA"))
sma(idat$TL_mm ~ idat$M_radius)
sma_plot <- idat |> 
  distinct(id, TL_mm, M_radius) |> 
  ggplot(aes(x = M_radius, y = TL_mm)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(intercept = standardized_major_axis[1, 1], slope = standardized_major_axis[2, 1],
              linewidth = 1, color = top_palette[1]) +
  labs(x = "Total otolith radius (mm)",
       y = "Total body length (mm)",
       title = "Total otolith radius vs. body length (Standardized Major Axis Regression)") +
  theme_minimal()

ggsave(sma_plot, file = here("figures", "sma_plot.png"))

# sma plot by location
get_sma_plot_points <- function(location_name) {
  data <- idat |> 
    filter(Location == location_name)
  sma <- with(data, line.cis(y = TL_mm, x = M_radius, method = "SMA"))
  intercept <- sma[1,1]
  slope <- sma[2,1]
  return(tibble("sma_intercept" = intercept,
                "sma_slope" = slope))
}

sma_plot_points <- purrr:::map_dfr(as.character(unique(idat$Location)), get_sma_plot_points) |> 
  mutate(Location = as.character(unique(idat$Location)))

sma_plot_by_loc <- idat |> 
  distinct(id, TL_mm, M_radius, Location) |>
  left_join(sma_plot_points, by = "Location") |> 
  ggplot(aes(x = M_radius, y = TL_mm)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(aes(intercept = sma_intercept, slope = sma_slope),
              linewidth = 1, color = top_palette[1]) +
  facet_wrap(~Location) +
  labs(x = "Total otolith radius (mm)",
       y = "Total body length (mm)",
       title = "Total otolith radius vs. body length (Standardized Major Axis Regression)") +
  theme_minimal()

ggsave(sma_plot_by_loc, file = here("figures", "sma_plot_by_loc.png"))



# plot otolith radius on total body length
otolith_to_body_length_plot <- idat |> 
  distinct(id, TL_mm, M_radius, Location) |> 
  ggplot(aes(x = M_radius, y = TL_mm)) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~Location) +
  labs(x = "Total otolith radius (mm)",
       y = "Total body length (mm)",
       title = "Total otolith radius vs. body length (Standardized Major Axis Regression)") +
  theme_minimal()

ggsave(otolith_to_body_length_plot, file = here("figures", "otolith_to_body_length_plot.png"))


# define back-calculation formulas ----------------------------------------

# biological intercept model Campana, pulled from Vigliola & Meekan 2009)
# bi_pars is a vector of length 2 with parameters for the model
calculate_biological_intercept <- function(capture_length, radius_increment, capture_radius, bi_pars) {
  length_increment <- capture_length + ((radius_increment - capture_radius) * ((capture_length - bi_pars[1]) / (capture_radius - bi_pars[2])))
  length_increment
}

# fl_pars is just b0
calculate_fraser_lee <- function(capture_length, fl_par, radius_increment, radius_capture) {
  length_increment <- fl_par + (capture_length - fl_par) * (radius_increment / radius_capture)
  length_increment
}

# fit bi_pars using sum of squares
bi_sum_of_squares <- function(bi_pars, data){
  data <- data |> 
    mutate(bi_length = calculate_biological_intercept(TL_mm, Radius, M_radius, bi_pars))
  
  sum_of_squares <- sum((data$bi_length - data$L)^2)
  return(sum_of_squares)
}

# minimize sum of squared differences
bi_starting_pars = c(L0 = 0, R0 = 0) # starting values
bi_fit <- nlminb(bi_starting_pars, bi_sum_of_squares, data = idat)
bi_fit$par
bi_fit$objective

L0p <- bi_fit$par[1]
R0p <- bi_fit$par[2]

bi_pars = c(L0 = L0p, R0 = R0p)
fl_par = 30 #b0, from MNDNR handbook
data_with_bc_lengths <- idat |> 
  mutate(bi_length = calculate_biological_intercept(TL_mm, Radius, M_radius, bi_pars),
         fl_length = calculate_fraser_lee(TL_mm, fl_par, Radius, M_radius))

data_with_bc_lengths |> 
  pivot_longer(c(bi_length, fl_length), names_to = "bc_method", values_to = "length_increment") |> 
  ggplot(aes(x = Age, y = length_increment, color = bc_method, 
             group = interaction(id, bc_method))) +
  geom_line(alpha = 0.2) +
  scale_color_manual(values = short_palette) +
  theme_minimal() +
  labs(x = "Age",
       y = "Back-calculated Lengths",
       title = "Fraser-Lee vs. Biological Intercept BCLs")

data_with_bc_lengths |> 
  ggplot(aes(x = Age, y = Radius, 
             group = id)) +
  geom_line(alpha = 0.4, linewidth = 0.9, color = short_palette[1]) +
  theme_minimal() +
  labs(x = "Age",
       y = "Otolith increment (mm)",
       title = "Individual otolith growth trajectories")

data_with_bc_lengths |> 
  mutate(residuals_bc_methods = bi_length - fl_length) |> 
  ggplot(aes(x = Age, y = residuals_bc_methods)) +
  geom_point(size = 2, alpha = 0.2, col =  short_palette[2]) +
  ylim(c(-60, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Age", 
       y = "BI BCL - FL BCL",
       title = "Back-Calculated Length (BCL) residuals for different methods") +
  theme_minimal()
  
data_with_bc_lengths |> 
  filter(Age == 1) |> 
  pivot_longer(c(bi_length, fl_length), names_to = "bc_method", values_to = "length_increment") |> 
  mutate("Back Calculation Method" = ifelse(bc_method == "bi_length", "Biological Intercept", "Fraser-Lee")) |> 
  ggplot(aes(x = length_increment , y = `Back Calculation Method`, fill = `Back Calculation Method`)) +
  geom_density_ridges(scale = 1, alpha = 0.8) +
  scale_fill_manual(values = short_palette) +
  labs(x = "Back-Calculated Length (mm)",
       title = "Back-Calculated Length Distribution for Age 1, Fraser-Lee vs. Biological Intercept") +
  theme_minimal()


# in previous script, added Fraser-Lee back-calculated lengths to idat as
# L_FL column and saved as 

save(idat, file='incdat_2021-BCM1.Rdata')


# plot back-calculated lengths vs. observed -------------------------------
binned_BCL <- idat |> 
  rename(age_at_capture = C_Age,
         age = Age,
         length_at_capture = TL_mm,
         BCL = L) |> 
  mutate(age_bin = case_when(age %in% 5:10 ~ "Ages 5-10",
                             age %in% 11:15 ~ "Ages 11-15",
                             age %in% 16:20 ~ "Ages 16-20",
                             age %in% 21:25 ~ "Ages 21-25",
                             age < 5 ~ "Ages 1-4",
                             TRUE ~ "Ages 26+")) |> 
  group_by(age_bin) |> 
  mutate(mean_BCL_length_by_bin = mean(BCL, na.rm = T),
         mean_length_at_capture = mean(length_at_capture, na.rm = T)) |> 
  pivot_longer(c(mean_BCL_length_by_bin, mean_length_at_capture),
               names_to = "mean_length_type", values_to = "mean_length") |> 
  glimpse()

age_bin_means <- binned_BCL |> 
  filter(mean_length_type == "mean_length_at_capture",
         age_bin != "Ages 1-4") |> 
  mutate(`Age Bin` = factor(age_bin, levels = c("Ages 5-10", "Ages 11-15", "Ages 16-20", "Ages 21-25", "Ages 26+"))) |> 
  distinct(`Age Bin`, mean_length)

age_bin_BCL_means <- binned_BCL |> 
  filter(mean_length_type == "mean_BCL_length_by_bin",
         age_bin != "Ages 1-4") |> 
  mutate(`Age Bin` = factor(age_bin, levels = c("Ages 5-10", "Ages 11-15", "Ages 16-20", "Ages 21-25", "Ages 26+"))) |> 
  distinct(`Age Bin`, mean_length)

binned_BCL |> 
  filter(age_bin != "Ages 1-4") |> 
  mutate(age_bin = factor(age_bin, levels = c("Ages 5-10", "Ages 11-15", "Ages 16-20", "Ages 21-25", "Ages 26+"))) |> 
  ggplot(aes(x = BCL, y = age_bin)) + 
  ggridges::geom_density_ridges(scale = 1, alpha = 0.5) +
  geom_vline(data = age_bin_means, 
             aes(xintercept = mean_length, color = age_bin),
             linewidth = 1) +
  scale_fill_manual(values = top_palette) +
  scale_color_manual(values = top_palette)


age_bin_BCL_plot <- binned_BCL |> 
  filter(age_bin != "Ages 1-4") |> 
  mutate(`Age Bin` = factor(age_bin, levels = c("Ages 5-10", "Ages 11-15", "Ages 16-20", "Ages 21-25", "Ages 26+"))) |> 
  ggplot(aes(x = BCL)) + 
  geom_density(aes(fill = `Age Bin`), alpha = 0.8) +
  geom_vline(data = age_bin_means, 
             aes(xintercept = mean_length), linetype = "dashed",
             linewidth = 1) +
  geom_vline(data = age_bin_BCL_means, 
             aes(xintercept = mean_length), linetype = "dotted",
             linewidth = 1) +
  facet_wrap(~`Age Bin`) +
  scale_fill_manual(values = top_palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Back-calculated length (mm)",
       y = "Density",
       title = "Back-calculated lengths by age bin")

ggsave(age_bin_BCL_plot, file = here("figures", "age_bin_BCL_plot.png"))

binned_BCL |> 
  ggplot(aes(x = age_bin, y = length, color = type)) + 
  geom_point()
