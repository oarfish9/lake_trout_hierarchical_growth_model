# plot residuals to show model fit

# plotting process error by year
library(tidyverse)
library(ggplot2)
library(here)

load(here("data", "best_fitting_rep.Rdata"))
load(here("data-raw", "incdat_2021.Rdata"))

top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")

short_palette = c("#354823", "#D69C4E")

pred_vs_observed <- tibble("Age" = rep$age,
                           "BCL_length" = rep$L,
                           "pred_length" = rep$L_hat,
                           "log_PE" = rep$log_PE,
                           "year_capture" = idat$year_capture,
                           "idx" = idat$id,
                           "pop_idx" = idat$pop_idx) |> 
  mutate(residuals = (BCL_length - pred_length),
         log_residuals = log(pred_length) - log(BCL_length),
         index = row_number())

pred_vs_observed |> 
  ggplot(aes(x = index, y = residuals)) +
  geom_point(color = top_palette[2], alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Row number",
       y = "Residuals",
       title = "Residuals of model-predicted lengths vs. observed (back-calculated) lengths")

pred_vs_observed |> 
  ggplot(aes(x = residuals)) +
  geom_density(fill = top_palette[2]) + 
  theme_minimal() +
  labs(x = "Residuals",
       y = "Density",
       title = "Residuals of model-predicted lengths vs. observed (back-calculated) lengths")

pred_vs_observed |> 
  ggplot(aes(x = BCL_length, y = pred_length)) +
  geom_point(color = top_palette[2], alpha = 0.2) +
  theme_minimal() +
  labs(x = "Back-calculated length",
       y = "Model-predicted length",
       title = "Relationship between back-calculated and predicted lengths")
summary(pred_vs_observed$residuals)


# year effects ------------------------------------------------------------

# plot observed (back-calculated) lengths vs. year-at-capture at age 4, age 15, and age 40 

BCL_by_year_age_plot <- pred_vs_observed |>  
  mutate(fake_date = as.Date(paste0(year_capture, "-01-01"))) |> 
  arrange(Age) |> 
  filter(Age %in% c(5, 10, 15, 20)) |> 
  mutate(Age = factor(paste0("Age ", Age), levels = c("Age 5", "Age 10", "Age 15", "Age 20"))) |> 
  ggplot(aes(x = fake_date, y = BCL_length)) +
  geom_point(color = top_palette[1], alpha = 0.7, size = 2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  facet_wrap(~Age) +
  theme_bw() +
  labs(x = "Year at capture", y = "Observed (back-calculated) length",
       title = "Back-calculated lengths at ages 5, 10, 20, and 25")

ggsave(BCL_by_year_age_plot, file = here("figures", "BCL_by_year_age_plot.png"))
