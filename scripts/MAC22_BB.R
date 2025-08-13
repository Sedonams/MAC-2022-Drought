
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

mycorrhizal_vars <- c(
  "amf_in_dry_soil", "rlc_p", "dse_in_dry_soil", "amf_tot", "dse_tot", "spores",
  "dse_hyphae", "fine_amf", "coarse_amf", "ves_and_arb", "dse_and_am",
  "vesicle_or_spore", "dse", "lse", "am_hyphae", "no_fungus", "arb","olpidium", "tot"
)

n_treatment <- length(unique(ds$treatment))
n_genotype <- length(unique(ds$genotype))
k_groups <- n_treatment * n_genotype

alpha <- 0.05
power <- 0.8

results_df <- data.frame(
  Variable = mycorrhizal_vars,
  f = NA_real_,
  n_per_group = NA_real_
)

for (i in seq_along(mycorrhizal_vars)) {
  var_name <- mycorrhizal_vars[i]
  cat("\nProcessing:", var_name, "\n")
  
  tmp <- ds[, c(var_name, "treatment", "genotype")]
  tmp <- tmp[complete.cases(tmp), ]
  
  tmp$group <- interaction(tmp$treatment, tmp$genotype)
  
  # Drop sparse groups
  group_counts <- table(tmp$group)
  cat("Group counts:", paste(names(group_counts), group_counts, collapse=", "), "\n")
  keep_groups <- names(group_counts[group_counts >= 2])
  tmp <- tmp[tmp$group %in% keep_groups, ]
  
  if (nrow(tmp) < 4 || length(unique(tmp$group)) < 2) {
    cat("Skipped: too few data after dropping sparse groups.\n")
    next
  }
  
  aov_model <- aov(tmp[[var_name]] ~ group, data = tmp)
  anova_tab <- anova(aov_model)
  
  if (!"group" %in% rownames(anova_tab)) {
    cat("Skipped: 'group' not in ANOVA table.\n")
    next
  }
  
  ss_effect <- anova_tab["group", "Sum Sq"]
  ss_resid  <- anova_tab["Residuals", "Sum Sq"]
  eta_sq <- ss_effect / (ss_effect + ss_resid)
  
  if (is.na(eta_sq) || eta_sq <= 0 || eta_sq >= 1) {
    cat("Skipped: invalid eta_sq.\n")
    next
  }
  
  f_val <- sqrt(eta_sq / (1 - eta_sq))
  k_val <- length(unique(tmp$group))
  
  pwr_result <- tryCatch(
    pwr.anova.test(
      k = k_val,
      f = f_val,
      sig.level = alpha,
      power = power
    ),
    error = function(e) {
      cat("Power calculation failed:", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(pwr_result)) {
    results_df$f[i] <- f_val
    results_df$n_per_group[i] <- ceiling(pwr_result$n)
  }
}

write.csv(results_df, "power_analysis_mycorrhizal.csv", row.names = FALSE)
cat("Power analysis results written to power_analysis_mycorrhizal.csv\n")
