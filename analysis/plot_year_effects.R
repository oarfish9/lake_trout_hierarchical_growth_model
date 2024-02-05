# plotting process error by year
library(tidyverse)
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
                                  "pop_idx" = idat$pop_idx,
                                  "PE" = exp(rep$log_PE))

random_fish_sample <- as.numeric(sample(process_error_with_year$id, 25))

process_error_with_year |> 
  filter(id %in% random_fish_sample) |> 
  ggplot(aes(x = year, y = PE)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Year",
       y = "Estimated Process Error",
       title = "Estimated Process Error by Year and Individual Fish") +
  facet_wrap(~id, scales = "free", ncol = 5) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme_minimal()



