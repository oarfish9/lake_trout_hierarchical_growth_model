# plotting process error by year
library(tidyverse)
library(here)
library(ggplot2)

load(here("data", "best_fitting_rep.Rdata"))
load(here("data-raw", "incdat_2021.Rdata"))

top_palette <- c("#1B0C42FF", "#FB9A06FF", "#781C6DFF","#CF4446FF", "#4B0C6BFF", "#ED6925FF", 
                 "#FCFFA4FF", "#000004FF", "#F7D03CFF", "#A52C60FF")

short_palette = c("#354823", "#D69C4E")


# extract process error estimates by individual and year ------------------

process_error_with_year <- tibble("log_PE" = rep$log_PE,
                                  "year" = idat$Year,
                                  "id" = as.numeric(idat$id),
                                  "pop_idx" = idat$pop_idx)

random_fish_sample <- as.numeric(sample(process_error_with_year$id, 25))

PE_with_year_plot <- process_error_with_year |> 
  filter(id %in% random_fish_sample) |> 
  mutate(fake_date = as.Date(paste0(year, "-01-01")),
         id = paste0("Individual ", id)) |> 
  ggplot(aes(x = fake_date, y = log_PE)) +
  geom_point(size = 2, alpha = 0.8, color = top_palette[1]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year",
       y = "Estimated log process error",
       title = "Estimated Process Error by Year and Individual Fish") +
  facet_wrap(~id, ncol = 5) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme_bw()

ggsave(PE_with_year_plot, file = here("figures", "PE_with_year_plot.png"))



