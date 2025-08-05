
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

#This is my path - change yours to whatever yours is
setwd("C:/Users/beabo/OneDrive/Documents/NAU/Sorghum/MAC-2022-Drought")


# Read and clean the data

#Where is sample ID in dmgr?
dmgr <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "DMGR Data")%>%
  clean_names()%>%
  rename(treatment = treat)

amf <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "AMF Score")%>%
  clean_names()%>%
  rename(treatment = treat)%>%
  mutate(genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
  position = as.integer(position))

mac <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "MAC22 Weights")%>%
  clean_names()%>%
  rename(treatment = treat)%>%
  filter(genotype != "No Tag")%>%
  mutate(rep = stringr::str_extract(sample_id, "[A-Za-z]"),
        position = stringr::str_extract(sample_id, "[0-9]+"),
        genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
        position = as.integer(position))

emh <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "EMH")%>%
  clean_names()%>%
  rename(sample_id = sample)%>%
  mutate(genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
         position = as.integer(position))%>%
  select(!geno) #Getting rid of geno column since I dont' know how it's different from the genotype column.

#Why does emh have both genotype and geno.

nuts <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "MAC22 Nutrients")%>%
  clean_names()%>%
  mutate(rep = stringr::str_extract(sample_id, "[A-Za-z]"),
        position = stringr::str_extract(sample_id, "[0-9]+"),
        genotype = stringr::str_replace(as.character(genotype), "\\.0$", ""),
        position = as.integer(position))

origins <- read_xlsx("data/MAC 2022 Mother for R.xlsx", sheet = "Geno Origins")%>%
  clean_names()

#What is G2C, and should I include it?
#g2c <- read_xlsx("data/G2C All Data (7).xlsx", sheet = "Weights")%>%
 # clean_names()
#g2camf <- read_xlsx("data/G2C All Data (7).xlsx", sheet = "AMF Score")%>%
#  clean_names()

#Same question for legacy study
#dsLegacy <- read_xlsx("C:/Users/sms689/Downloads/Legacy Study 2024.xlsx", sheet = "LP2")

#What are dsdrought and dswet?
#dsdrought<- read_xlsx("C:/Users/sms689/Downloads/z6-19231(z6-19231)-1736291922.xlsx", sheet = "Config 1")
#dswet<- read_xlsx("C:/Users/sms689/Downloads/z6-19230(z6-19230)-1736199516.xlsx", sheet = "Config 1")



#-----------Combining Data----

ds <- mac %>%
  left_join(amf, by = c("genotype", "treatment", "rep","position", "sample_id")) %>%
  left_join(emh, by = c("genotype", "treatment", "rep","position", "sample_id")) %>%
  left_join(nuts,by = c("genotype", "treatment", "rep","position", "sample_id"))%>%
  mutate(across(c(florets, shoot_wt, main_shoot_wtkg, amf_in_dry_soil), as.numeric))

#what about dmgr? leaving out for now - does not have sample ID

sapply(ds, class)

#Skipping g2c again for now bc I don't know what it is
#g2c$Rep <- as.integer(g2c$Rep)
#g2c$Genotype <- gsub("\\.0$", "", as.character(g2c$Genotype))
#g2camf$Genotype <- gsub("\\.0$", "", as.character(g2camf$Genotype))

#g2combo <- g2c %>%  
 # left_join(g2camf, by = c("Genotype", "Treat", "Rep"))


write.csv(ds, "data/MAC22_cleaned.csv", row.names = FALSE)

ds %>% group_by(genotype)%>%tally()%>%filter(n<5) #Cool.
