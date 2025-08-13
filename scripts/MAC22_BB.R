
#BB
#8/5/25
#This script was originally written by SS. I'm adapting it for some analyses for the Sorghum symposium happening in 2 weeks.
#Doing some visualizations.

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
library(forcats)


#Set up working directory
getwd()

#This is my path - change yours to whatever yours is
setwd("C:/Users/beabo/OneDrive/Documents/NAU/Sorghum/MAC-2022-Drought")

ds <- read.csv("data/MAC22_cleaned.csv") #Using the cleaned up dataset. Missing some things like dmgr. 

theme_set(theme_bw())

#------macombo plots----

vars_og <- c("shoot_wt", "florets","arb","vesicle_or_spore", "amf_in_dry_soil" , "rlc_p","dse_in_dry_soil", "d13c", "c", "d15n", "n", "p", "cn_ratio", "np_ratio")

vars <- c("arb", "am_hyphae", "dse", "lse", "vesicle_or_spore", "amf_in_dry_soil" , "rlc_p","dse_in_dry_soil")

ds_scaled <- ds %>%
  group_by(genotype) %>%
  mutate(across(all_of(vars), ~ scale(.)[,1])) %>%
  ungroup()

# Step 2: Safe t-tests for each genotype × variable
stats_df <- ds_scaled %>%
  select(genotype, treatment, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by(genotype, variable) %>%
  summarise(
    p_value = tryCatch({
      if (n_distinct(treatment) == 2) {
        t.test(value ~ treatment)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    stars = case_when(
      !is.na(p_value) & p_value <= 0.001 ~ "***",
      !is.na(p_value) & p_value <= 0.01  ~ "**",
      !is.na(p_value) & p_value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Step 3: Calculate drought effect (for coloring)
heatmap_df <- ds_scaled %>%
  select(genotype, treatment, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by(genotype, treatment, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = mean_val) %>%
  mutate(effect = Droughted - Watered)

# Step 4: Order variables from most red to most blue
var_order <- heatmap_df %>%
  group_by(variable) %>%
  summarise(mean_effect = mean(effect, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_effect)) %>%
  pull(variable)

heatmap_df <- heatmap_df %>%
  mutate(variable = factor(variable, levels = var_order))

# Step 5: Merge stats into heatmap data
heatmap_df <- left_join(heatmap_df, stats_df, by = c("genotype", "variable"))

# Step 6: Plot heatmap with stars
print(
  ggplot(heatmap_df, aes(x = variable, y = genotype, fill = effect)) +
    geom_tile(color = "white") +
    geom_text(aes(label = stars), color = "black", size = 4) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Std. diff\n(Drought - Watered)"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      x = "Variable",
      y = "Genotype",
      title = "Drought effects by genotype (standardized)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      plot.title = element_text(hjust = 0.5)
    )
)


#Power analysis: how many samples need to be scored




# Power analysis for mycorrhizal data
library(pwr)

# Identify relevant mycorrhizal variables (replace with actual variable names)
mycorrhizal_vars <- c("amf_in_dry_soil", "rlc_p", "dse_in_dry_soil") #Example - replace with your actual variables

# Create a data frame to store results
results_df <- data.frame(Variable = character(), n = numeric(), stringsAsFactors = FALSE)

# Perform power analysis (example using t-test; adjust for ANOVA if needed)
for (var in mycorrhizal_vars) {
  # Assuming two treatment groups (adjust if more)
  power_result <- pwr.t.test(
    d = 0.5, # Effect size (Cohen's d) - adjust based on prior knowledge or pilot data
    sig.level = 0.05,
    power = 0.8, # Desired power (80% is common)
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Variable = var, n = power_result$n))
}

# Write results to CSV
tryCatch({
  write.csv(results_df, "power_analysis_results.csv", row.names = FALSE)
  cat("Power analysis results written to power_analysis_results.csv\n")
}, error = function(e) {
  cat(paste0("Error writing CSV: ", e$message), "\n")
})

# Interpretation:
# n: The required sample size per group to achieve the specified power.
# If n is larger than your current sample size, you need to collect more samples.


