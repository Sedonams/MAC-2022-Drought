
#BB
#8/5/25
#This script was originally written by SS. I'm adapting it for some analyses for the Sorghum symposium happening in 2 weeks.
#Doing some visualizations.
#editing Bea's charts

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
library(ggfortify)  # For autoplot
library(scales)
library(viridis)
library(ggrepel)  # For smart label positioning


#Set up working directory
getwd()

#This is my path - change yours to whatever yours is
setwd("H:/MAC-2022-Drought")

ds <- read.csv("data/MAC22_cleaned.csv") #Using the cleaned up dataset. Missing some things like dmgr. 

theme_set(theme_bw())

#------Drought Response Heatmap for nutrients biomass and florets----

vars <- c("shoot_wt", "florets", "d13c", "c", "d15n", "n", "p", "cn_ratio", "np_ratio")

# Step 1: Scale data
ds_scaled <- ds %>%
  group_by(genotype) %>%
  mutate(across(all_of(vars), ~ scale(.)[,1])) %>%
  ungroup()

# Step 2: Create long format data
filtered_long_data <- ds_scaled %>%
  select(genotype, treatment, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))  # Remove NAs

# Step 3: Safe t-tests for all genotype-variable combinations
stats_df <- filtered_long_data %>%
  group_by(genotype, variable) %>%
  summarise(
    p_value = tryCatch({
      if (length(unique(treatment)) == 2) {
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

# Step 4: Calculate drought effect for all combinations
heatmap_df <- filtered_long_data %>%
  group_by(genotype, treatment, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = mean_val) %>%
  mutate(effect = Droughted - Watered) %>%
  filter(!is.na(effect))  # Remove any remaining NAs

# Step 5: Order variables from most red to most blue
var_order <- heatmap_df %>%
  group_by(variable) %>%
  summarise(mean_effect = mean(effect, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_effect)) %>%
  pull(variable)

heatmap_df <- heatmap_df %>%
  mutate(variable = factor(variable, levels = var_order))

# Step 6: Merge stats into heatmap data
heatmap_df <- left_join(heatmap_df, stats_df, by = c("genotype", "variable"))

# Step 7: Enhanced Plot
library(scales)   # For better formatting

# Create prettier variable names for display
pretty_names <- c(
  "shoot_wt" = "Shoot Weight",
  "florets" = "Florets",
  "arb" = "Arbuscules",
  "vesicle_or_spore" = "Vesicles/Spores",
  "amf_in_dry_soil" = "AMF in Dry Soil",
  "amf_rlc" = "Root Length Colonized",
  "dse_in_dry_soil" = "DSE in Dry Soil",
  "d13c" = "δ¹³C",
  "c" = "Carbon Content",
  "d15n" = "δ¹⁵N",
  "n" = "Nitrogen Content",
  "p" = "Phosphorus Content",
  "cn_ratio" = "C:N Ratio",
  "np_ratio" = "N:P Ratio"
)

# Apply prettier names
heatmap_df <- heatmap_df %>%
  mutate(variable_pretty = factor(pretty_names[as.character(variable)], 
                                  levels = pretty_names[var_order]))

# Create the enhanced plot
p <- ggplot(heatmap_df, aes(x = variable_pretty, y = genotype, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = stars), color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", 
    midpoint = 0,
    name = "Standardized\nDifference\n(Drought - Watered)",
    breaks = pretty_breaks(n = 5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL,
    y = "Genotype",
    title = "Drought Response Across Genotypes",
    subtitle = "Standardized differences between drought and watered treatments",
    caption = "Significance: *** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05"
  ) +
  theme_minimal() +
  theme(
    # Text styling
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
    
    # Axis styling
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    
    # Legend styling
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(1.2, "cm"),
    
    # Panel styling
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # Plot margins
    plot.margin = margin(20, 20, 20, 20),
    
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p)

# Optional: Save high-quality version
ggsave("drought_response_heatmap_simple.png", plot = p, width = 14, height = 8, dpi = 300, bg = "white")





#------Making RLC's for the colonization data----

ds <- ds %>%
  mutate(arb_ves_and_arb = rowSums(select(., arb, ves_and_arb), na.rm = TRUE)) %>%
  mutate(vesicle_or_spore_ves_and_arb = rowSums(select(., vesicle_or_spore, ves_and_arb), na.rm = TRUE)) %>%
  mutate(am_hyphae_dse_and_am = rowSums(select(., am_hyphae, dse_and_am), na.rm = TRUE)) %>%
  mutate(dse_dse_and_am = rowSums(select(., dse, dse_and_am), na.rm = TRUE))

# Count non-zero values in each column to confirm successful combination
sum(ds$arb_ves_and_arb != 0, na.rm = TRUE)
sum(ds$arb != 0, na.rm = TRUE)

columns_to_convert <- c("am_hyphae_dse_and_am","arb_ves_and_arb", "dse_dse_and_am", "vesicle_or_spore_ves_and_arb", "lse", "coil",
"olpidium",	"mold", "plasmodiophorids",	"dot_line",	"non_am",	"fine_endo")

ds <- ds %>%
  mutate(across(all_of(columns_to_convert), 
                ~ if_else(!is.na(.x) & !is.na(tot), 
                          (.x / tot) * 100, 
                          NA_real_),
                .names = "{.col}_rlc"))

ds <- ds %>%
  rename(amf_rlc = rlc_p)

#------Drought Response Heatmap for Colonization----

vars <- c("amf_rlc", "am_hyphae_dse_and_am_rlc","arb_ves_and_arb_rlc", "dse_dse_and_am_rlc",
          "vesicle_or_spore_ves_and_arb_rlc", "lse_rlc", "coil_rlc",
          "olpidium_rlc",	"mold_rlc", "plasmodiophorids_rlc",	"dot_line_rlc",	"non_am_rlc",	
          "fine_endo_rlc","amf_in_dry_soil","dse_in_dry_soil","shoot_wt", "florets", "d13c",
          "c", "d15n", "n", "p", "cn_ratio", "np_ratio")

# Step 1: Scale data
ds_scaled <- ds %>%
  group_by(genotype) %>%
  mutate(across(all_of(vars), ~ scale(.)[,1])) %>%
  ungroup()

# Step 2: Identify valid genotypes based on amf_rlc column (need 3+ reps per treatment)
valid_genotypes <- ds_scaled %>%
  filter(!is.na(amf_rlc)) %>%  # Only look at non-NA values in amf_rlc
  group_by(genotype, treatment) %>%
  summarise(n_reps = n(), .groups = "drop") %>%  # Count reps per treatment
  pivot_wider(names_from = treatment, values_from = n_reps, values_fill = 0) %>%
  filter(Droughted >= 3 & Watered >= 3) %>%  # Both treatments must have 3+ reps
  pull(genotype)

print(paste("Keeping genotypes:", paste(valid_genotypes, collapse = ", ")))

# Step 3: Filter all data to only valid genotypes
filtered_long_data <- ds_scaled %>%
  filter(genotype %in% valid_genotypes) %>%
  select(genotype, treatment, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))  # Remove NAs after filtering

# Step 4: Safe t-tests for valid genotypes only
stats_df <- filtered_long_data %>%
  group_by(genotype, variable) %>%
  summarise(
    p_value = tryCatch({
      if (length(unique(treatment)) == 2) {
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

# Step 5: Calculate drought effect for valid genotypes only
heatmap_df <- filtered_long_data %>%
  group_by(genotype, treatment, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = mean_val) %>%
  mutate(effect = Droughted - Watered) %>%
  filter(!is.na(effect))  # Remove any remaining NAs

# Step 6: Order variables from most red to most blue
var_order <- heatmap_df %>%
  group_by(variable) %>%
  summarise(mean_effect = mean(effect, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_effect)) %>%
  pull(variable)

heatmap_df <- heatmap_df %>%
  mutate(variable = factor(variable, levels = var_order))

# Step 7: Merge stats into heatmap data
heatmap_df <- left_join(heatmap_df, stats_df, by = c("genotype", "variable"))

# Step 8: Enhanced Plot
library(viridis)  # For better color palettes
library(scales)   # For better formatting

# Create prettier variable names for display
pretty_names <- c(
  "amf_rlc" = "AMF RLC",
  "am_hyphae_dse_and_am_rlc" = "AM Hyphae RLC",
  "arb_ves_and_arb_rlc" = "Arbuscules RLC",
  "dse_dse_and_am_rlc" = "DSE RLC",
  "vesicle_or_spore_ves_and_arb_rlc" = "Vesicles/Spores RLC",
  "lse_rlc" = "LSE RLC",
  "coil_rlc" = "Coils RLC",
  "olpidium_rlc" = "Olpidium RLC",
  "mold_rlc" = "Mold RLC",
  "plasmodiophorids_rlc" = "Plasmodiophorids RLC",
  "dot_line_rlc" = "Dot Line RLC",
  "non_am_rlc" = "Non-AM RLC",
  "fine_endo_rlc" = "Fine Endophytes RLC",
  "amf_in_dry_soil" = "AMF in Dry Soil",
  "dse_in_dry_soil" = "DSE in Dry Soil",
  "shoot_wt" = "Shoot Weight",
  "florets" = "Florets",
  "arb" = "Arbuscules",
  "vesicle_or_spore" = "Vesicles/Spores",
  "amf_in_dry_soil" = "AMF in Dry Soil",
  "amf_rlc" = "Root Length Colonized",
  "dse_in_dry_soil" = "DSE in Dry Soil",
  "d13c" = "δ¹³C",
  "c" = "Carbon Content",
  "d15n" = "δ¹⁵N",
  "n" = "Nitrogen Content",
  "p" = "Phosphorus Content",
  "cn_ratio" = "C:N Ratio",
  "np_ratio" = "N:P Ratio"
)

# Apply prettier names
heatmap_df <- heatmap_df %>%
  mutate(variable_pretty = factor(pretty_names[as.character(variable)], 
                                  levels = pretty_names[var_order]))

# Create the enhanced plot
p <- ggplot(heatmap_df, aes(x = variable_pretty, y = genotype, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = stars), color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", 
    midpoint = 0,
    name = "Standardized\nDifference\n(Drought - Watered)",
    breaks = pretty_breaks(n = 5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL,
    y = "Genotype",
    title = "Drought Response Across Genotypes",
    subtitle = "Standardized differences between drought and watered treatments",
    caption = "Significance: *** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05\nOnly genotypes with ≥3 replicates per treatment shown"
  ) +
  theme_minimal() +
  theme(
    # Text styling
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
    
    # Axis styling
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    
    # Legend styling
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(1.2, "cm"),
    
    # Panel styling
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # Plot margins
    plot.margin = margin(20, 20, 20, 20),
    
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p)

# Optional: Save high-quality version
ggsave("drought_response_heatmap.png", plot = p, width = 14, height = 8, dpi = 300, bg = "white")





#Power analysis: how many samples need to be scored----



# Power analysis for mycorrhizal data
library(pwr)

mycorrhizal_vars <- c(
  "amf_in_dry_soil", "amf_rlc", "dse_in_dry_soil", "amf_tot", "dse_tot", "spores",
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


#---Ordination (PCA) for plant traits and leaf nutrients----
# Define all variables for Ordination, plant traits and leaf nutrients
all_vars <- c("shoot_wt", "florets", "d13c",
  "c", "d15n", "n", "p")

# Prepare data for PCA
pca_data <- ds %>%
  select(genotype, treatment, all_of(all_vars)) %>%
  # Remove rows with too many missing values (optional threshold)
  filter(rowSums(is.na(select(., all_of(all_vars)))) <= length(all_vars) * 0.5) %>%
  # Create a unique identifier for each sample
  mutate(sample_id = paste(genotype, treatment, row_number(), sep = "_"))

# Remove variables that are mostly NA (optional)
vars_to_keep <- pca_data %>%
  select(all_of(all_vars)) %>%
  summarise_all(~ sum(!is.na(.)) / n()) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "prop_complete") %>%
  filter(prop_complete >= 0.3) %>%  # Keep variables with at least 30% complete data
  pull(variable)

print(paste("Keeping", length(vars_to_keep), "variables for PCA"))
print(paste("Variables kept:", paste(vars_to_keep, collapse = ", ")))

# Final data preparation Just for Treatment
#pca_data_final <- pca_data %>%
 # select(sample_id, genotype, treatment, all_of(vars_to_keep)) %>%
  # Remove any remaining rows with all NA values
 # filter(!if_all(all_of(vars_to_keep), is.na))


# Final data preparation with Genotype X Treatment
pca_data_final <- pca_data %>%
  select(genotype, treatment, all_of(vars_to_keep)) %>%
  group_by(genotype, treatment) %>%
  summarise(across(all_of(vars_to_keep), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  mutate(sample_id = paste(genotype, treatment, sep = "_")) %>%
  filter(!if_all(all_of(vars_to_keep), is.na))


# Perform PCA (handling missing values)
pca_matrix <- pca_data_final %>%
  select(all_of(vars_to_keep)) %>%
  as.matrix()

# Option 1: Use only complete cases (most conservative)
complete_cases <- complete.cases(pca_matrix)
if(sum(complete_cases) < 10) {
  cat("Too few complete cases, using imputation method...\n")
  
  # Option 2: Simple mean imputation for missing values
  for(i in 1:ncol(pca_matrix)) {
    col_mean <- mean(pca_matrix[, i], na.rm = TRUE)
    pca_matrix[is.na(pca_matrix[, i]), i] <- col_mean
  }
  
  pca_matrix_clean <- pca_matrix
  pca_data_clean <- pca_data_final
  
} else {
  # Use complete cases
  cat(paste("Using", sum(complete_cases), "complete cases out of", nrow(pca_matrix), "total samples\n"))
  pca_matrix_clean <- pca_matrix[complete_cases, ]
  pca_data_clean <- pca_data_final[complete_cases, ]
}

# Check for any remaining issues
if(any(is.na(pca_matrix_clean)) || any(is.infinite(pca_matrix_clean))) {
  # Additional cleaning if needed
  pca_matrix_clean[is.na(pca_matrix_clean)] <- 0
  pca_matrix_clean[is.infinite(pca_matrix_clean)] <- 0
}

# Perform PCA using prcomp
pca_result <- prcomp(pca_matrix_clean, scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
pc3_var <- round(var_explained[3], 1)

# Create PCA data frame for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  genotype = pca_data_clean$genotype,
  treatment = pca_data_clean$treatment,
  sample_id = pca_data_clean$sample_id
)

# Get loadings for arrows
loadings_df <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[, 1],
  PC2 = pca_result$rotation[, 2],
  PC3 = pca_result$rotation[, 3]
) %>%
  # Scale loadings for better visualization
  mutate(
    PC1_scaled = PC1 * 3,
    PC2_scaled = PC2 * 3,
    PC3_scaled = PC3 * 3
  )

# Create prettier variable names for loadings
pretty_var_names <- c(
  "amf_rlc" = "AMF RLC",
  "am_hyphae_dse_and_am_rlc" = "AM Hyphae",
  "arb_ves_and_arb_rlc" = "Arbuscules",
  "dse_dse_and_am_rlc" = "DSE",
  "vesicle_or_spore_ves_and_arb_rlc" = "Vesicles/Spores",
  "lse_rlc" = "LSE",
  "coil_rlc" = "Coils",
  "olpidium_rlc" = "Olpidium",
  "mold_rlc" = "Mold",
  "plasmodiophorids_rlc" = "Plasmodiophorids",
  "dot_line_rlc" = "Dot Line",
  "non_am_rlc" = "Non-AM",
  "fine_endo_rlc" = "Fine Endophytes",
  "amf_in_dry_soil" = "AMF in Dry Soil",
  "dse_in_dry_soil" = "DSE in Dry Soil",
  "shoot_wt" = "Shoot Weight",
  "florets" = "Florets",
  "d13c" = "δ¹³C",
  "c" = "Carbon Content",
  "d15n" = "δ¹⁵N",
  "n" = "Nitrogen Content",
  "p" = "Phosphorus Content",
  "cn_ratio" = "C:N Ratio",
  "np_ratio" = "N:P Ratio"
)

loadings_df <- loadings_df %>%
  mutate(Variable_pretty = ifelse(Variable %in% names(pretty_var_names), 
                                  pretty_var_names[Variable], 
                                  Variable))

# Create the main PCA plot, Selective labeling with variable arrow thickness 
#Only label the most important loadings and use different strategies for others
# Calculate loading strength for arrow thickness
loadings_with_strength <- loadings_df %>% 
  filter(abs(PC1) > 0.3 | abs(PC2) > 0.3) %>%
  mutate(
    # Calculate the magnitude of loading vector
    loading_strength = sqrt(PC1^2 + PC2^2),
    # Create thickness categories
    thickness_category = case_when(
      loading_strength >= 0.7 ~ "strong",      # Thickest arrows
      loading_strength >= 0.5 ~ "medium",     # Medium arrows
      TRUE ~ "weak"                           # Thinnest arrows
    ),
    # Map to actual linewidth values
    arrow_thickness = case_when(
      thickness_category == "strong" ~ 1.2,
      thickness_category == "medium" ~ 0.8,
      thickness_category == "weak" ~ 0.4
    ),
    # Arrow alpha based on strength
    arrow_alpha = case_when(
      thickness_category == "strong" ~ 0.9,
      thickness_category == "medium" ~ 0.7,
      thickness_category == "weak" ~ 0.5
    )
  )

#code for the plot
p1_selective <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = treatment), 
               alpha = 0.3, type = "norm", level = 0.68, linewidth = 1) +
  geom_point(aes(color = treatment, shape = treatment), size = 1.5, alpha = 0.8) +
  # Arrows with variable thickness based on loading strength
  geom_segment(data = loadings_with_strength,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled,
                   linewidth = arrow_thickness, alpha = arrow_alpha),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "gray20") +
  # Custom linewidth scale
  scale_linewidth_identity() +
  scale_alpha_identity() +
  # Label only the strongest loadings directly
  geom_text_repel(data = loadings_with_strength %>% 
                    filter(loading_strength > 0.5), # Use loading_strength instead
                  aes(x = PC1_scaled * 1.1, y = PC2_scaled * 1.1, 
                      label = Variable_pretty),
                  size = 3.5, 
                  color = "black", 
                  fontface = "bold",
                  box.padding = 0.6,
                  point.padding = 0.3,
                  force = 3) +
  # Add full labels for medium-strength loadings using ggrepel
  geom_text_repel(data = loadings_with_strength %>% 
                    filter(loading_strength <= 0.5 & loading_strength > 0.3),
                  aes(x = PC1_scaled * 1.05, y = PC2_scaled * 1.05, 
                      label = Variable_pretty), # Full labels, not abbreviated
                  size = 2.8, 
                  color = "gray40", 
                  fontface = "plain",
                  box.padding = 0.3,
                  point.padding = 0.2,
                  force = 1.5,
                  max.overlaps = Inf) +
  scale_color_manual(values = c("Droughted" = "#D73027", "Watered" = "#1A9850"),
                     name = "Treatment") +
  scale_shape_manual(values = c("Droughted" = 17, "Watered" = 16),
                     name = "Treatment") +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "Principal Component Analysis (Treatment X Genotype)",
    subtitle = "Plant performance and leaf chemistry",
    caption = "Arrow thickness indicates loading strength"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.caption = element_text(size = 9, color = "gray50")
  )

# creat plot
print(p1_selective)

# Print loading strength summary for reference
cat("\nLoading strength summary:\n")
strength_summary <- loadings_with_strength %>%
  select(Variable_pretty, loading_strength, thickness_category) %>%
  arrange(desc(loading_strength))
print(strength_summary)


# Create a scree plot
scree_data <- data.frame(
  PC = paste0("PC", 1:min(10, ncol(pca_result$x))),
  Variance = var_explained[1:min(10, length(var_explained))]
) %>%
  mutate(PC = factor(PC, levels = PC))

p2 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3, fontface = "bold") +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "Scree Plot",
    subtitle = "Variance explained by each PC"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p2)


# Print summary statistics
cat("\n=== PCA SUMMARY ===\n")
cat(paste("Total samples included:", nrow(pca_df), "\n"))
cat(paste("Total variables included:", length(vars_to_keep), "\n"))
cat(paste("PC1 variance explained:", pc1_var, "%\n"))
cat(paste("PC2 variance explained:", pc2_var, "%\n"))
cat(paste("PC3 variance explained:", pc3_var, "%\n"))
cat(paste("Cumulative variance (PC1+PC2):", round(pc1_var + pc2_var, 1), "%\n"))

# Optional: Save plots
ggsave("pca_biplot.png", plot = p1_selective, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("pca_scree.png", plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

#
#
#
#
#




#---Ordination (PCA) for colonization (as well as plant traits and leaf nutrients)----
# Colonization vars, plant traits, and leaf nutrients
all_vars <- c(
  "amf_rlc", "am_hyphae_dse_and_am_rlc", "arb_ves_and_arb_rlc", "dse_dse_and_am_rlc",
  "vesicle_or_spore_ves_and_arb_rlc", "lse_rlc", "coil_rlc",
  "olpidium_rlc", "mold_rlc", "plasmodiophorids_rlc", "dot_line_rlc", "non_am_rlc",
  "fine_endo_rlc", "amf_in_dry_soil", "dse_in_dry_soil", "shoot_wt", "florets", "d13c",
  "c", "d15n", "n", "p")

# Prepare data for PCA
pca_data <- ds %>%
  select(genotype, treatment, all_of(all_vars)) %>%
  # Remove rows with too many missing values (optional threshold)
  filter(rowSums(is.na(select(., all_of(all_vars)))) <= length(all_vars) * 0.5) %>%
  # Create a unique identifier for each sample
  mutate(sample_id = paste(genotype, treatment, row_number(), sep = "_"))

# Remove variables that are mostly NA (optional)
vars_to_keep <- pca_data %>%
  select(all_of(all_vars)) %>%
  summarise_all(~ sum(!is.na(.)) / n()) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "prop_complete") %>%
  filter(prop_complete >= 0.3) %>%  # Keep variables with at least 30% complete data
  pull(variable)

print(paste("Keeping", length(vars_to_keep), "variables for PCA"))
print(paste("Variables kept:", paste(vars_to_keep, collapse = ", ")))

# Final data preparation Just for Treatment
#pca_data_final <- pca_data %>%
# select(sample_id, genotype, treatment, all_of(vars_to_keep)) %>%
# Remove any remaining rows with all NA values
# filter(!if_all(all_of(vars_to_keep), is.na))


# Final data preparation with Genotype X Treatment
pca_data_final <- pca_data %>%
  select(genotype, treatment, all_of(vars_to_keep)) %>%
  group_by(genotype, treatment) %>%
  summarise(across(all_of(vars_to_keep), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  mutate(sample_id = paste(genotype, treatment, sep = "_")) %>%
  filter(!if_all(all_of(vars_to_keep), is.na))


# Perform PCA (handling missing values)
pca_matrix <- pca_data_final %>%
  select(all_of(vars_to_keep)) %>%
  as.matrix()

# Option 1: Use only complete cases (most conservative)
complete_cases <- complete.cases(pca_matrix)
if(sum(complete_cases) < 10) {
  cat("Too few complete cases, using imputation method...\n")
  
  # Option 2: Simple mean imputation for missing values
  for(i in 1:ncol(pca_matrix)) {
    col_mean <- mean(pca_matrix[, i], na.rm = TRUE)
    pca_matrix[is.na(pca_matrix[, i]), i] <- col_mean
  }
  
  pca_matrix_clean <- pca_matrix
  pca_data_clean <- pca_data_final
  
} else {
  # Use complete cases
  cat(paste("Using", sum(complete_cases), "complete cases out of", nrow(pca_matrix), "total samples\n"))
  pca_matrix_clean <- pca_matrix[complete_cases, ]
  pca_data_clean <- pca_data_final[complete_cases, ]
}

# Check for any remaining issues
if(any(is.na(pca_matrix_clean)) || any(is.infinite(pca_matrix_clean))) {
  # Additional cleaning if needed
  pca_matrix_clean[is.na(pca_matrix_clean)] <- 0
  pca_matrix_clean[is.infinite(pca_matrix_clean)] <- 0
}

# Perform PCA using prcomp
pca_result <- prcomp(pca_matrix_clean, scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
pc3_var <- round(var_explained[3], 1)

# Create PCA data frame for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  genotype = pca_data_clean$genotype,
  treatment = pca_data_clean$treatment,
  sample_id = pca_data_clean$sample_id
)

# Get loadings for arrows
loadings_df <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[, 1],
  PC2 = pca_result$rotation[, 2],
  PC3 = pca_result$rotation[, 3]
) %>%
  # Scale loadings for better visualization
  mutate(
    PC1_scaled = PC1 * 3,
    PC2_scaled = PC2 * 3,
    PC3_scaled = PC3 * 3
  )

# Create prettier variable names for loadings
pretty_var_names <- c(
  "amf_rlc" = "AMF RLC",
  "am_hyphae_dse_and_am_rlc" = "AM Hyphae",
  "arb_ves_and_arb_rlc" = "Arbuscules",
  "dse_dse_and_am_rlc" = "DSE",
  "vesicle_or_spore_ves_and_arb_rlc" = "Vesicles/Spores",
  "lse_rlc" = "LSE",
  "coil_rlc" = "Coils",
  "olpidium_rlc" = "Olpidium",
  "mold_rlc" = "Mold",
  "plasmodiophorids_rlc" = "Plasmodiophorids",
  "dot_line_rlc" = "Dot Line",
  "non_am_rlc" = "Non-AM",
  "fine_endo_rlc" = "Fine Endophytes",
  "amf_in_dry_soil" = "AMF in Dry Soil",
  "dse_in_dry_soil" = "DSE in Dry Soil",
  "shoot_wt" = "Shoot Weight",
  "florets" = "Florets",
  "d13c" = "δ¹³C",
  "c" = "Carbon Content",
  "d15n" = "δ¹⁵N",
  "n" = "Nitrogen Content",
  "p" = "Phosphorus Content",
  "cn_ratio" = "C:N Ratio",
  "np_ratio" = "N:P Ratio"
)

loadings_df <- loadings_df %>%
  mutate(Variable_pretty = ifelse(Variable %in% names(pretty_var_names), 
                                  pretty_var_names[Variable], 
                                  Variable))

# Create Selective labeling with variable arrow thickness 
#Only label the most important loadings and use different strategies for others

# Calculate loading strength for arrow thickness
loadings_with_strength <- loadings_df %>% 
  filter(abs(PC1) > 0.3 | abs(PC2) > 0.3) %>%
  mutate(
    # Calculate the magnitude of loading vector
    loading_strength = sqrt(PC1^2 + PC2^2),
    # Create thickness categories
    thickness_category = case_when(
      loading_strength >= 0.7 ~ "strong",      # Thickest arrows
      loading_strength >= 0.5 ~ "medium",     # Medium arrows
      TRUE ~ "weak"                           # Thinnest arrows
    ),
    # Map to actual linewidth values
    arrow_thickness = case_when(
      thickness_category == "strong" ~ 1.2,
      thickness_category == "medium" ~ 0.8,
      thickness_category == "weak" ~ 0.4
    ),
    # Arrow alpha based on strength
    arrow_alpha = case_when(
      thickness_category == "strong" ~ 0.9,
      thickness_category == "medium" ~ 0.7,
      thickness_category == "weak" ~ 0.5
    )
  )



genotypes_to_plot <- c("156203", "157033", "181080", "181083", "641815", "E29W1")

pca_df_filtered <- pca_df %>%
  filter(genotype %in% genotypes_to_plot)

# Compute ellipse label positions (center of each treatment group)
ellipse_labels <- pca_df_filtered %>%
  group_by(treatment) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")






p1_selective <- ggplot(pca_df_filtered, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = treatment), 
               alpha = 0.3, type = "norm", level = 0.68, linewidth = 1) +
  geom_point(aes(color = genotype, shape = treatment), size = 3, alpha = 0.8)+
  
  
  geom_segment(data = loadings_with_strength,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled,
                   linewidth = arrow_thickness, alpha = arrow_alpha),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "gray20") +
  # Custom linewidth scale
  scale_linewidth_identity() +
  scale_alpha_identity() +
  # Label only the strongest loadings directly
  geom_text_repel(data = loadings_with_strength %>% 
                    filter(loading_strength > 0.5), # Use loading_strength instead
                  aes(x = PC1_scaled * 1.1, y = PC2_scaled * 1.1, 
                      label = Variable_pretty),
                  size = 3.5, 
                  color = "black", 
                  fontface = "bold",
                  box.padding = 0.6,
                  point.padding = 0.3,
                  force = 3) +
  # Add full labels for medium-strength loadings using ggrepel
  geom_text_repel(data = loadings_with_strength %>% 
                    filter(loading_strength <= 0.5 & loading_strength > 0.3),
                  aes(x = PC1_scaled * 1.05, y = PC2_scaled * 1.05, 
                      label = Variable_pretty), # Full labels, not abbreviated
                  size = 2.8, 
                  color = "gray40", 
                  fontface = "plain",
                  box.padding = 0.3,
                  point.padding = 0.2,
                  force = 1.5,
                  max.overlaps = Inf) +
                  scale_color_discrete(name = "Genotype") +
                  #scale_color_manual(values = c(
                   # "Droughted" = "#D73027",
                  #  "Watered" = "darkblue" ))
                  scale_shape_manual(values = c("Droughted" = 17, "Watered" = 16),
                   name = "Treatment") +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "Principal Component Analysis (Preliminary AMF Colonization Genotypes)",
    subtitle = "Plant performance, chemistry, and fungal colonization traits",
    caption = "Arrow thickness indicates loading strength"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.caption = element_text(size = 9, color = "gray50")
  )

# Display the selective approach with variable arrow thickness
print(p1_selective)








# Print loading strength summary for reference
cat("\nLoading strength summary:\n")
strength_summary <- loadings_with_strength %>%
  select(Variable_pretty, loading_strength, thickness_category) %>%
  arrange(desc(loading_strength))
print(strength_summary)


# Create a scree plot
scree_data <- data.frame(
  PC = paste0("PC", 1:min(10, ncol(pca_result$x))),
  Variance = var_explained[1:min(10, length(var_explained))]
) %>%
  mutate(PC = factor(PC, levels = PC))

p2 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3, fontface = "bold") +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "Scree Plot",
    subtitle = "Variance explained by each PC"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p2)


# Print summary statistics
cat("\n=== PCA SUMMARY ===\n")
cat(paste("Total samples included:", nrow(pca_df), "\n"))
cat(paste("Total variables included:", length(vars_to_keep), "\n"))
cat(paste("PC1 variance explained:", pc1_var, "%\n"))
cat(paste("PC2 variance explained:", pc2_var, "%\n"))
cat(paste("PC3 variance explained:", pc3_var, "%\n"))
cat(paste("Cumulative variance (PC1+PC2):", round(pc1_var + pc2_var, 1), "%\n"))

# Optional: Save plots
ggsave("pca_biplot.png", plot = p1_selective, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("pca_scree.png", plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

genotypes_to_plot <- c("156203", "157033", "181080", "181083", "641815", "E29W1")

pca_df_filtered <- pca_df %>%
  filter(genotype %in% genotypes_to_plot)

# Compute ellipse label positions (center of each treatment group)
ellipse_labels <- pca_df_filtered %>%
  group_by(treatment) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

p1_treatment_ellipse <- ggplot(pca_df_filtered, aes(x = PC1, y = PC2)) +
  # Ellipses colored by treatment
  stat_ellipse(aes(color = treatment), 
               alpha = 0.3, type = "norm", level = 0.68, linewidth = 1.2) +
  
  # Label treatment ellipses directly
  geom_text_repel(data = ellipse_labels,
                  aes(x = PC1, y = PC2, label = treatment, color = treatment),
                  size = 4.5, fontface = "bold",
                  segment.color = NA) +
  
  # Points: colored by genotype, shaped by treatment
  geom_point(aes(color = genotype, shape = treatment), size = 2.5, alpha = 0.9) +
  
  # Arrows (loadings)
  geom_segment(data = loadings_with_strength,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled,
                   linewidth = arrow_thickness, alpha = arrow_alpha),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray20") +
  
  # Labels for strong loadings
  geom_text_repel(data = loadings_with_strength %>% filter(loading_strength > 0.5),
                  aes(x = PC1_scaled * 1.1, y = PC2_scaled * 1.1, label = Variable_pretty),
                  size = 3.5, color = "black", fontface = "bold",
                  box.padding = 0.6, point.padding = 0.3, force = 3) +
  
  # Manual color and shape scales
  scale_color_manual(values = c(
    "Droughted" = "#D73027",
    "Watered" = "#1A9850"
  )) +
  
  scale_shape_manual(values = c("Droughted" = 17, "Watered" = 16)) +
  
  # Remove redundant legends
  guides(color = "none") +
  
  # Labels and theme
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA: Ellipses by Treatment, Points Colored by Genotype",
    subtitle = "Colonization, Traits & Nutrients | Genotypes: selected subset",
    caption = "Ellipses show 68% confidence by treatment; arrows = loading strength"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.caption = element_text(size = 9, color = "gray50")
  )
# creat plot
print(p1_treatment_ellipse)
