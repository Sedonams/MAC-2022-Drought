
#BB
#8/5/25
#This script was originally written by SS. I'm adapting it for some analyses for the Sorghum symposium happening in 2 weeks.
#Cleaning up the main dataset for analysis, so it's in good shape and everyone's using the same version.

#-----Loading Libraries and Data Sets----
# Load necessary libraries
library(readxl)   # For reading Excel files
library(janitor)  # For cleaning column names
library(dplyr)    # For using the pipe operator (%>%)
library(corrr)    # For nicer correlation matrices
library(tidyr)
library(ggplot2)
library(EnvStats)
library(lubridate)
library(stringr)
library(gridExtra)


#Set up working directory
getwd()
#This is my path
setwd("H:/MAC-2022-Drought")

# Read and clean the data
dmgr <- read_xlsx("data/MAC 2022 Mother.xlsx", sheet = "DMGR Data") %>%
  clean_names() %>% # nolint
  rename(treatment = treat) %>%
  rename(age = days_from_sprout_or_transplant_to_harvest) %>%
  mutate(age = as.integer(age))
# convert dmgr data to date safely
date_cols <- c("seed_or_reseed_date", "sprout_date",
               "transplant_date", "harvest_date")
for (col in date_cols) {
  if (col %in% names(dmgr)) {
    # Specify the date format, e.g. "%Y-%m-%d" or "%m/%d/%Y" as needed
    dmgr[[col]] <- as.Date(as.character(dmgr[[col]]),
                           format = "%Y-%m-%dT%H:%M:%OSZ")
  }
}

colnames(dmgr)

# onto the amf scoring tab cleanup
columns_to_convert <- c("am_hyphae", "arb", "dse", "vesicle", "lse",
                        "coil", "fine_endo", "olpidium", "mold",
                        "plasmodiophorids", "dot_line", "non_am")

amf <- read_xlsx("data/MAC 2022 Mother.xlsx", sheet = "AMF Score") %>%
  clean_names() %>%
  rename(treatment = treat) %>%
  mutate(genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
         position = as.integer(position),
         non_am = rowSums(dplyr::select(., olpidium, mold,
                                        plasmodiophorids, dot_line),
                          na.rm = TRUE),
         arb = rowSums(dplyr::select(., arb, ves_and_arb),
                       na.rm = TRUE),
         vesicle =
           rowSums(dplyr::select(., vesicle_or_spore, ves_and_arb),
                   na.rm = TRUE),
         am_hyphae =
           rowSums(dplyr::select(., course_am_hyphae, dse_and_am),
                   na.rm = TRUE),
         dse = rowSums(dplyr::select(., dse, dse_and_am),
                       na.rm = TRUE)) %>%
  mutate(across(ends_with("_rlc"), as.numeric)) %>%  # Ensure numeric
  mutate(across(all_of(columns_to_convert),
                ~ if_else(!is.na(.x) & !is.na(tot),
                          (.x / tot) * 100,
                          NA_real_),
                .names = "{.col}_rlc")) %>%
  mutate(across(all_of(columns_to_convert), # Replace zero values with NA
                ~ na_if(.x, 0))) %>%
  mutate(across(ends_with("_rlc"), #Replace NA with 0 where tot > 0
                ~ if_else(tot > 0 & is.na(.x), 0, .x)))


colnames(amf)
View(amf)

#now cleaning up the main weights data
mac <- read_xlsx("data/MAC 2022 Mother.xlsx", sheet = "MAC22 Weights") %>%
  clean_names() %>%
  rename(treatment = treat) %>%
  filter(genotype != "No Tag") %>%
  mutate(rep = stringr::str_extract(sample_id, "[A-Za-z]"),
         position = stringr::str_extract(sample_id, "[0-9]+"),
         genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
         position = as.integer(position))

colnames(mac)


#cleaning up the Extramatical Hyphae (emh) data
emh <- read_xlsx("data/MAC 2022 Mother.xlsx", sheet = "EMH") %>%
  clean_names() %>%
  mutate(genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
         position = as.integer(position))

colnames(emh)

#Leaf Tissue Nutrient Data
nuts <- read_xlsx("data/MAC 2022 Mother.xlsx",
                  sheet = "MAC22 Nutrients") %>%
  clean_names() %>%
  mutate(rep = stringr::str_extract(sample_id, "[A-Za-z]"),
         position = stringr::str_extract(sample_id, "[0-9]+"),
         genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
         position = as.integer(position))


#Genotype Origin data
origins <- read_xlsx("data/MAC 2022 Mother.xlsx",
                     sheet = "Geno Origins") %>%
  clean_names()

#Moisture probe data
soilmoisture <- read_xlsx("data/mac22_all_probes.xlsx") %>%
  clean_names() %>%
  mutate(timestamps = as.POSIXct(timestamps, format = "%Y-%m-%dT%H:%M:%OSZ",
                                 tz = "UTC"))

colnames(soilmoisture)
View(soilmoisture)

#-----------Combining Data----
ds <- mac  %>%
  left_join(amf, by = c("genotype", "treatment", "rep",
                        "position", "sample_id")) %>%
  left_join(emh, by = c("genotype", "treatment", "rep",
                        "position", "sample_id")) %>%
  left_join(nuts, by = c("genotype", "treatment", "rep",
                         "position", "sample_id")) %>%
  mutate(across(c(florets, shoot_wt, main_shoot_wtkg, amf_in_dry_soil, tot),
                as.numeric))

glimpse(ds)
View(ds)

write.csv(ds, "data/MAC22_cleaned.csv", row.names = FALSE)

ds %>% group_by(genotype) %>% tally() %>% filter(n < 5) #Cool. what is this
