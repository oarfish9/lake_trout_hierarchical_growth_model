library(readxl)
library(tidyverse)
library(here)
library(janitor)

# spatial analyses
survey_metadata <- read_xlsx("data-raw/SuperiorSurveyData.xlsx",
                             sheet = "Surveys")

surveys_with_centroids <- survey_metadata |> 
  clean_names() |>
  distinct(location, year, latitude, longitude)


load("data-raw/incdat_2021.Rdata")

fish_ids_with_centroids <- idat |> 
  left_join(surveys_with_centroids, 
            by = c("Location" = "location", 
                   "year_capture" = "year")) |> 
  distinct(id, latitude, longitude, location = Location, year_capture)

save(fish_ids_with_centroids, file = "data/fish_ids_with_centroids.Rdata")

