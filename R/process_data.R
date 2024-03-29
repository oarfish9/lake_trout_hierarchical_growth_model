# read in back-calculation data (recieved from Mike Hansen) to include only
# Lean morphs
# create indexing variable and year of capture based on serial number
# refactor of data_wrangling_05-11-2021

# load libraries
library(tidyverse)
library(here)
library(readxl)

# read in increment data
increment_data_raw <- read.table(here("data-raw", "mydata_2021.txt"), 
                                 header = T, sep = '\t')
increment_data <- increment_data_raw |> 
  filter(Morph == "Lean")

# read in data with biological data
biodat_raw <- read_xlsx(here("data-raw", "Lake Trout Bio Data.xlsx"))

# process to only the individuals we're interested in
biodat <- biodat_raw |> 
  select(Maturity, Sex, id = ID, Air_g, Wet_g)

# now join
full_data_raw <- merge(biodat, increment_data, by = "id")

# create indexing column
full_data <- full_data_raw |> 
  mutate(id = factor(id),
         Location = factor(Location),
         Serial = sub("_", "-", Serial),
         ycode = substr(Serial, 1, regexpr("\\-", Serial) - 1),
         year_capture = case_when(ycode == "LSBR06" ~ 2006,
                                  ycode %in% c("LSBR14", "LSSR14") ~ 2014,
                                  ycode == "LS02" ~ 2002,
                                  ycode == "LS03" ~ 2003,
                                  ycode == "LSIR06" ~ 2006,
                                  ycode == "LSIR07" ~ 2007,
                                  ycode == "LS04" ~ 2004,
                                  ycode %in% c("LSSR13", "LSSS13") ~ 2013)) |> 
  group_by(id) |> 
  mutate(idx1 = cur_group_id(), # indexing for R
         idx = idx1- 1, # change indexing to start at 0 for C++
         max_age = max(Age),
         year_factor = row_number() - 1) |> 
  ungroup() |> 
  mutate(Year = year_capture - (max_age - year_factor) + 1,
         pop_idx = cur_group_id(), .by = Location,
         pop_idx = pop_idx - 1) |> 
  ungroup() |> 
  select(-c(max_age, year_factor))

idat <- full_data
save(idat, file = here("data-raw", "incdat_2021.Rdata"))

biodat_save <- biodat
save(biodat_save, file = here("data-raw", "incdat_2022_biodat2.Rdata"))

# locations
id_col <- full_data |> 
  distinct(Location, pop_idx)

write.csv(id_col, here("data", "Lake Trout Lake Names and IDs.csv"))

# documentation of raw data
colnames(idat)
idat$TL_mm # length at capture, nind unique values
idat$M_radius # otolith radius at capture, nind unique values
idat$Radius # radius at each age, calculated by Mike using BI model
