# back-calculation methods

# TODO refactor

rm(list=ls())
# libraries
library(dplyr)
library(ggplot2)
library(here)
library(wesanderson)


# plot otolith v body -----------------------------------------------------
load(here("data-raw", "incdat_2021.Rdata"))
short_palette = c("#354823", "#D69C4E")

# we don't have the individual trajectories; we can only show
# fit of our linear model to the data

idat |> 
  distinct(id, TL_mm, M_radius) |> 
  ggplot(aes(x = M_radius, y = TL_mm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Total otolith radius (mm)",
       y = "Total body length (mm)",
       title = "Total otolith radius vs. body length") +
  theme_minimal()

total_length_radius_linear <- lm(TL_mm ~ M_radius, data = idat)


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