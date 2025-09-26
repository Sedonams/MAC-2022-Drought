# Mixed Linear Models for dmgr and MAC22

#install.packages(c("readr","dplyr","stringr","lme4","lmerTest","broom.mixed","emmeans","performance","ggplot2"))
library(readr)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(performance)
library(ggplot2)

getwd()

#This is my path -
setwd("H:/MAC-2022-Drought")

MAC22_CLEANED <- read.csv("data/MAC22_cleaned.csv") #Using the cleaned up dataset. Missing some things like dmgr. 
MAC22_DMGR <- read.csv("data/mac22_dmgr.csv") #Using the cleaned up dataset. Missing some things like dmgr. 
MAC22_ORIGINS <- read.csv("data/mac22_origins.csv") #Using the cleaned up dataset. Missing some things like dmgr. 
MAC22_OUT <- read.csv("data") #Using the cleaned up dataset. Missing some things like dmgr. 



#Response = biomass (shoot, root, leaf), nutrient content, or fungal colonization.
#Fixed effects = treatment, genotype origin traits, genotype
#Random effects = NumofPlants (since multiple samples per genotype), rep, position, or scorer.

MAC22_MLM_ShootWt <- ShootWt ~ Treat * genotype + (1|NumofPlants) + (1|Rep)





