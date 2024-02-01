library(readxl)
library(tidyverse)
library(here)
library(janitor)
library(ggplot2)

# spatial analyses
survey_metadata <- read_xlsx("data-raw/SuperiorSurveyData.xlsx",
                             sheet = "Surveys")

surveys_with_centroids <- survey_metadata |> 
  clean_names() |>
  distinct(location, year, latitude, longitude) |> 
  mutate(longitude = -1 * longitude)


load("data-raw/incdat_2021.Rdata")

fish_ids_with_centroids <- idat |> 
  left_join(surveys_with_centroids, 
            by = c("Location" = "location", 
                   "year_capture" = "year")) |> 
  distinct(id, latitude, longitude, location = Location, year_capture)

save(fish_ids_with_centroids, file = "data/fish_ids_with_centroids.Rdata")


# get individual net locations
net_metadata <- read_xlsx("data-raw/SuperiorSurveyData.xlsx",
                          sheet = "Sampling") |> 
  clean_names()

nets_with_lat_longs <- net_metadata |> 
  mutate(longitude = -1 * longitude) |> 
  distinct(location, survey_id, effort_id, latitude, longitude) |> 
  glimpse()

save(nets_with_lat_longs, file = "data/nets_with_lat_longs.Rdata")

# read in outline of lake superior

# Lake Superior: https://www.glerl.noaa.gov/data/bathy/bathy.html
lake_superior_lat_longs <- read_csv("data-raw/LakeSuperior_latlong.csv") |> 
  rename(longitude = `long...1`, latitude = `long...2`, location = Location) |> 
  glimpse()

save(lake_superior_lat_longs, file = "data/lake_superior_lat_longs.Rdata")


