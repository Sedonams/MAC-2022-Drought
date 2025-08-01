#-----Loading Libraries and Data Sets----
# Load necessary libraries
library(readxl)   # For reading Excel files
library(janitor)  # For cleaning column names
library(dplyr)    # For using the pipe operator (%>%)
library(corrr)    # For nicer correlation matrices

# Read and clean the data
dmgr <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "DMGR Data")

amf <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "AMF Score")
mac <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "MAC22 Weights")
emh <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "EMH")
nuts <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "MAC22 Nutrients")
Origins <- read_xlsx("C:/Users/sms689/Downloads/MAC 2022 Mother (16).xlsx", sheet = "Geno Origins")

g2c <- read_xlsx("C:/Users/sms689/Downloads/G2C All Data (7).xlsx", sheet = "Weights")
g2camf <- read_xlsx("C:/Users/sms689/Downloads/G2C All Data (7).xlsx", sheet = "AMF Score")

dsLegacy <- read_xlsx("C:/Users/sms689/Downloads/Legacy Study 2024.xlsx", sheet = "LP2")

dsdrought<- read_xlsx("C:/Users/sms689/Downloads/z6-19231(z6-19231)-1736291922.xlsx", sheet = "Config 1")
dswet<- read_xlsx("C:/Users/sms689/Downloads/z6-19230(z6-19230)-1736199516.xlsx", sheet = "Config 1")


# Load unnecessary libraries
library(tidyr)
library(ggplot2)
library(EnvStats)
library(lubridate)
library(gridExtra)

#-----------Combining Data----

mac$Rep <- gsub("[^A-Za-z]", "", mac$`Sample ID`)  # Remove everything that's not a letter
mac$Position <- gsub("[^0-9]", "", mac$`Sample ID`)  # Remove everything that's not a digit

nuts$Rep <- gsub("[^A-Za-z]", "", nuts$`Sample ID`)  # Remove everything that's not a letter
nuts$Position <- gsub("[^0-9]", "", nuts$`Sample ID`)  # Remove everything that's not a digit

colnames(amf)[colnames(amf) == "Treatment"] <- "Treat" 
colnames(nuts)[colnames(nuts) == "Treatment"] <- "Treat"
colnames(emh)[colnames(emh) == "Treatment"] <- "Treat"

colnames(amf)[colnames(amf) == "Sample_ID"] <- "Sample"
colnames(nuts)[colnames(nuts) == "Sample ID"] <- "Sample" 
colnames(emh)[colnames(emh) == "Sample ID"] <- "Sample"
colnames(mac)[colnames(mac) == "Sample ID"] <- "Sample"

mac$Genotype <- gsub("\\.0$", "", as.character(mac$Genotype))
amf$Genotype <- gsub("\\.0$", "", as.character(amf$Genotype))
emh$Genotype <- gsub("\\.0$", "", as.character(emh$Genotype))
nuts$Genotype <- gsub("\\.0$", "", as.character(nuts$Genotype))

mac$Position <- as.integer(mac$Position)
nuts$Position <- as.integer(nuts$Position)
emh$Position <- as.integer(emh$Position)
amf$Position <- as.integer(amf$Position)

macombo <- mac %>%
  left_join(amf, by = c("Genotype", "Treat", "Rep","Position", "Sample")) %>%
  left_join(emh, by = c("Genotype", "Treat", "Rep","Position", "Sample")) %>%
  left_join(nuts,by = c("Genotype", "Treat", "Rep","Position", "Sample"))

macombo$Florets <- as.numeric(macombo$Florets)
macombo$ShootWt <- as.numeric(macombo$ShootWt)
macombo$MainShootWtkg <- as.numeric(macombo$MainShootWtkg)
macombo$AMFInDrySoil <- as.numeric(as.character(macombo$AMFInDrySoil))

g2c$Rep <- as.integer(g2c$Rep)
g2c$Genotype <- gsub("\\.0$", "", as.character(g2c$Genotype))
g2camf$Genotype <- gsub("\\.0$", "", as.character(g2camf$Genotype))

g2combo <- g2c %>%  
  left_join(g2camf, by = c("Genotype", "Treat", "Rep"))

#------macombo plots----

# Filter out 0 values and NA values
#subset the data in macombo into just the genotypes scored for emh
macombosubset <- macombo %>%
  filter(!is.na(TOT) & 
           TOT != 0 & 
           TOT != "#DIV/0!")

# Fit linear models for each treatment
model_watered <- lm(AMFInDrySoil ~ RLC_P, data = macombosubset %>% filter(Treat == "Watered"))
model_droughted <- lm(AMFInDrySoil ~ RLC_P, data = macombosubset %>% filter(Treat == "Droughted"))

# Extract slope and R-squared values
slope_watered <- coef(model_watered)[2]  # Slope for "Watered"
r2_watered <- summary(model_watered)$r.squared  # R-squared for "Watered"

slope_droughted <- coef(model_droughted)[2]  # Slope for "Droughted"
r2_droughted <- summary(model_droughted)$r.squared  # R-squared for "Droughted"





ds <- macombo %>%
  group_by(Genotype) %>%
  mutate(rr = if_else(Treat == "Droughted", 
                      log(ShootWt / mean(ShootWt[Treat == "Watered"], na.rm = TRUE)), 
                      NA_real_))

ds %>%
  group_by(Genotype) %>%
  mutate(mean_rr = mean(rr, na.rm = TRUE)) %>%  # Calculate mean rr for each Genotype
  ungroup() %>%
  arrange(desc(mean_rr)) %>%
  ggplot(aes(x = fct_reorder(factor(Genotype), rr, .fun = mean, .desc = TRUE), 
             y = rr, 
             fill = Genotype)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept = 0) +
  labs(title  = "Shoot Weight Response Ratio by Genotype", x = "Genotype", y = "Shoot Weight Response Ratio")+
  theme(legend.position = "none")

# Plot with faceting by NumOfPlants from the dmgr dataset
ggplot(data = ds,
       aes(x = Genotype, y = mgr)) +
  geom_boxplot(aes(fill = Genotype)) +  # Color boxplot
  theme(axis.title.x = element_blank()) +  # Remove X title
  facet_grid(~NumOfPlants)+ # Facet by NumOfPlants from the dmgr dataset
  stat_n_text() 




# Plot
macombosubset %>%
  ggplot(aes(x = RLC_P, y = AMFInDrySoil, fill = Genotype, color = Treat)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  scale_shape_manual(values = c("Watered" = 16, "Droughted" = 17)) +
  scale_color_manual(values = c("Watered" = "lightblue", "Droughted" = "darkred")) +
  geom_text(
    aes(x = 0.1, y = max(AMFInDrySoil) * 0.8, label = paste("Watered: Slope =", round(slope_watered, 2),
                                                            "\nR² =", round(r2_watered, 2))),
    color = "lightblue", hjust = 0, size = 3
  ) +
  geom_text(
    aes(x = 0.1, y = max(AMFInDrySoil) * 0.9, label = paste("Droughted: Slope =", round(slope_droughted, 2),
                                                            "\nR² =", round(r2_droughted, 2))),
    color = "darkred", hjust = 0, size = 3
  )







# Create the AMF boxplot# Create the AMF boxplotGenotype
macombosubset %>%
  ggplot(aes(x = Genotype, y = AMFInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "AMF in Dry Soil (m/g)")

# Create the DSE boxplot
macombosubset %>%
  ggplot(aes(x = Genotype, y = DSEInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "DSE in Dry Soil (m/g)")

#violin of Treat
macombosubset %>%
  ggplot(aes(x = Treat, y = AMFInDrySoil, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())



#------correlations with Origins data-----
# Correlation Matrices with Origins Data - Updated for Actual Column Names
# Load required libraries
library(dplyr)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(reshape2)
library(corrr)
library(RColorBrewer)

# 1. Prepare and Clean Origins Data

# Clean the Origins data
Origins$Genotype <- gsub("\\.0$", "", as.character(Origins$Genotype))

# Print available columns to confirm
print("=== ORIGINS DATASET COLUMNS ===")
print(colnames(Origins))

# Define continuous variables from Origins for correlation analysis
origins_continuous <- c(
  "Origin Anthesis Date (DAP)*",  # Note: may need to handle the asterisk
  "Anthesis date (days)",
  "Harvest date (days)", 
  "Total fresh weight (kg)",
  "Brix (maturity)",
  "Brix (milk)",
  "Dry weight (kg)", 
  "Stalk height (cm)",
  "Dry tons per acre",
  "ADF (% DM)",      # Acid Detergent Fiber
  "NDF (% DM)",      # Neutral Detergent Fiber  
  "NFC (% DM)",      # Non-Fiber Carbohydrates
  "Lignin (% DM)"
)

# Define categorical variables for grouping analyses
origins_categorical <- c(
  "Common Name",
  "Taxonomy", 
  "GRIN Origin",
  "Use",
  "Type",
  "Photoperiod", 
  "Subspecies",
  "Growth",
  "Taxa",
  "Pericarp",
  "pigmentation"
)

# Clean and convert continuous variables to numeric
origins_clean <- Origins
for(col in origins_continuous) {
  if(col %in% colnames(origins_clean)) {
    origins_clean[[col]] <- as.numeric(as.character(origins_clean[[col]]))
  }
}

# Combine with macombo data
macombo_origins <- macombo %>%
  left_join(origins_clean, by = "Genotype")

# 2. Plant Traits + Origins Continuous Variables Correlation

# Define plant trait variables
trait_vars <- c("ShootWt", "Florets", "d15N", "d13C", "N", "P", "RLC_P", 
                "AMFInDrySoil", "DSEInDrySoil")

# Get available continuous origins variables (remove those with all NAs)
available_origins_continuous <- origins_continuous[
  sapply(origins_continuous, function(x) {
    if(x %in% colnames(macombo_origins)) {
      sum(!is.na(macombo_origins[[x]])) > 5  # At least 5 non-NA values
    } else {
      FALSE
    }
  })
]

print("=== AVAILABLE CONTINUOUS ORIGINS VARIABLES ===")
print(available_origins_continuous)

# 3. Genotype-Level Correlation Matrix

# Calculate genotype means for traits (averaging across treatments and reps)
genotype_means <- macombo_origins %>%
  group_by(Genotype) %>%
  summarise(
    across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)),
    across(all_of(available_origins_continuous), ~ first(.x[!is.na(.x)])),
    .groups = 'drop'
  )

# Select complete cases for correlation
all_vars <- c(trait_vars, available_origins_continuous)
cor_data <- genotype_means %>%
  select(all_of(all_vars)) %>%
  select_if(function(x) sum(!is.na(x)) > 5) %>%  # Keep variables with >5 observations
  drop_na()

print(paste("Number of genotypes with complete data:", nrow(cor_data)))
print(paste("Variables included in correlation:", ncol(cor_data)))

# Function to create correlation heatmap with significance
create_detailed_heatmap <- function(data, title = "Correlation Matrix") {
  cor_result <- rcorr(as.matrix(data))
  
  cor_df <- melt(cor_result$r)
  pval_df <- melt(cor_result$P)
  names(cor_df) <- c("Var1", "Var2", "Correlation")
  names(pval_df) <- c("Var1", "Var2", "p_value")
  
  cor_full <- left_join(cor_df, pval_df, by = c("Var1", "Var2")) %>%
    mutate(
      sig_label = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**", 
        p_value <= 0.05 ~ "*",
        TRUE ~ ""
      ),
      label_text = ifelse(is.na(p_value), 
                          round(Correlation, 2),
                          paste0(round(Correlation, 2), sig_label)),
      # Create categories for better visualization
      var1_type = ifelse(Var1 %in% trait_vars, "Plant Trait", "Origin Trait"),
      var2_type = ifelse(Var2 %in% trait_vars, "Plant Trait", "Origin Trait")
    )
  
  ggplot(cor_full, aes(x = Var2, y = Var1, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = label_text), size = 2.5) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
    theme_minimal() +
    coord_fixed() +
    labs(title = title,
         subtitle = "Significance: *** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05",
         x = "", y = "") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 12),
      legend.position = "right"
    )
}

# Create main correlation heatmap
if(ncol(cor_data) > 1) {
  create_detailed_heatmap(cor_data, "Plant Traits + Origins Characteristics")
}

# 4. Focus on Plant Traits vs Origins Only

# Create subset focusing on cross-correlations between plant traits and origins
if(length(available_origins_continuous) > 0 && length(trait_vars) > 0) {
  
  # Calculate correlation matrix
  full_cor_matrix <- cor(cor_data, use = "complete.obs")
  
  # Extract plant traits vs origins correlations
  trait_origin_cors <- full_cor_matrix[trait_vars, available_origins_continuous, drop = FALSE]
  
  # Melt for plotting
  trait_origin_df <- melt(trait_origin_cors)
  names(trait_origin_df) <- c("Plant_Trait", "Origin_Trait", "Correlation")
  
  # Create focused heatmap
  ggplot(trait_origin_df, aes(x = Origin_Trait, y = Plant_Trait, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
    theme_minimal() +
    coord_fixed() +
    labs(title = "Plant Traits vs Origin Characteristics",
         x = "Origin Characteristics", y = "Plant Traits") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 5. Treatment-Specific Correlations with Origins

# Watered treatment correlations
watered_data <- macombo_origins %>%
  filter(Treat == "Watered") %>%
  group_by(Genotype) %>%
  summarise(
    across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)),
    across(all_of(available_origins_continuous), ~ first(.x[!is.na(.x)])),
    .groups = 'drop'
  ) %>%
  select(all_of(all_vars)) %>%
  select_if(function(x) sum(!is.na(x)) > 5) %>%
  drop_na()

# Droughted treatment correlations  
droughted_data <- macombo_origins %>%
  filter(Treat == "Droughted") %>%
  group_by(Genotype) %>%
  summarise(
    across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)),
    across(all_of(available_origins_continuous), ~ first(.x[!is.na(.x)])),
    .groups = 'drop'
  ) %>%
  select(all_of(all_vars)) %>%
  select_if(function(x) sum(!is.na(x)) > 5) %>%
  drop_na()

# Compare correlation matrices
if(nrow(watered_data) > 5 && nrow(droughted_data) > 5) {
  watered_cor <- cor(watered_data, use = "complete.obs")
  droughted_cor <- cor(droughted_data, use = "complete.obs")
  
  # Create side-by-side corrplots
  par(mfrow = c(1, 2))
  corrplot(watered_cor, method = "color", type = "upper", 
           title = "Watered Treatment", mar = c(0,0,2,0), tl.cex = 0.6)
  corrplot(droughted_cor, method = "color", type = "upper",
           title = "Droughted Treatment", mar = c(0,0,2,0), tl.cex = 0.6)
  par(mfrow = c(1, 1))
}

# 6. Response Ratios vs Origins 
# Calculate drought response ratios
response_ratios <- macombo_origins %>%
  group_by(Genotype) %>%
  summarise(
    across(all_of(trait_vars), 
           ~ log(mean(.x[Treat == "Droughted"], na.rm = TRUE) / 
                   mean(.x[Treat == "Watered"], na.rm = TRUE)),
           .names = "rr_{.col}"),
    across(all_of(available_origins_continuous), ~ first(.x[!is.na(.x)])),
    .groups = 'drop'
  )

# Response ratio variable names
rr_vars <- paste0("rr_", trait_vars)
available_rr_vars <- rr_vars[rr_vars %in% colnames(response_ratios)]

# Correlation of response ratios with origins
rr_data <- response_ratios %>%
  select(all_of(c(available_rr_vars, available_origins_continuous))) %>%
  select_if(function(x) sum(is.finite(x)) > 5) %>%
  filter(if_all(everything(), is.finite))

if(nrow(rr_data) > 5 && ncol(rr_data) > 1) {
  create_detailed_heatmap(rr_data, "Drought Response Ratios vs Origins")
}

# 7. Categorical Variable Analysis
# Analyze key categorical variables
key_categorical <- c("Type", "Use", "Subspecies", "Growth")
available_categorical <- key_categorical[key_categorical %in% colnames(macombo_origins)]

for(cat_var in available_categorical) {
  if(length(unique(macombo_origins[[cat_var]])) > 1 && 
     length(unique(macombo_origins[[cat_var]])) < 10) {
    
    print(paste("\n=== ANALYSIS BY", cat_var, "==="))
    
    # Group means by categorical variable
    group_means <- macombo_origins %>%
      filter(!is.na(.data[[cat_var]])) %>%
      group_by(.data[[cat_var]]) %>%
      summarise(
        n = n(),
        across(all_of(trait_vars), ~ mean(.x, na.rm = TRUE)),
        .groups = 'drop'
      )
    
    print(group_means)
  }
}

# 8. Summary and Key Findings

print("\n=== CORRELATION ANALYSIS SUMMARY ===")
print(paste("Total genotypes analyzed:", length(unique(macombo_origins$Genotype))))
print(paste("Plant traits included:", length(trait_vars)))
print(paste("Origins characteristics included:", length(available_origins_continuous)))

# Find strongest correlations between plant traits and origins
if(exists("trait_origin_cors") && length(trait_origin_cors) > 0) {
  print("\n=== STRONGEST TRAIT-ORIGIN CORRELATIONS ===")
  
  # Find top positive and negative correlations
  max_indices <- which(abs(trait_origin_cors) == max(abs(trait_origin_cors), na.rm = TRUE), arr.ind = TRUE)
  if(nrow(max_indices) > 0) {
    for(i in 1:min(5, nrow(max_indices))) {
      row_idx <- max_indices[i, 1]
      col_idx <- max_indices[i, 2]
      trait_name <- rownames(trait_origin_cors)[row_idx]
      origin_name <- colnames(trait_origin_cors)[col_idx]
      cor_value <- trait_origin_cors[row_idx, col_idx]
      print(paste(trait_name, "vs", origin_name, ":", round(cor_value, 3)))
    }
  }
  
  # Summary statistics
  print(paste("\nMean absolute correlation (traits vs origins):", 
              round(mean(abs(trait_origin_cors), na.rm = TRUE), 3)))
  print(paste("Max absolute correlation:", 
              round(max(abs(trait_origin_cors), na.rm = TRUE), 3)))
}

#------correlation matrix of Drought Response Ratios for ShootWt, Florets, d15N, d13C, N, P, and RLC_P########
library(dplyr)
library(tidyr)
library(Hmisc)
library(reshape2)
library(ggplot2)

# Define traits
traits <- c("ShootWt", "Florets", "d15N", "d13C", "N", "P", "RLC_P")

# Compute log response ratio
rr_data <- macombo %>%
  group_by(Genotype) %>%
  mutate(across(all_of(traits), 
                ~ if_else(Treat == "Droughted",
                          log(.x / mean(.x[Treat == "Watered"], na.rm = TRUE)),
                          NA_real_),
                .names = "rr_{.col}")) %>%
  ungroup()

# Extract RR variables and clean
rr_matrix <- rr_data %>%
  select(starts_with("rr_")) %>%
  filter(if_all(everything(), is.finite))  # Remove rows with NA, NaN, or Inf

# Compute correlations and p-values
rr_mat <- as.matrix(rr_matrix)
cor_res <- rcorr(rr_mat)

cor_df <- melt(cor_res$r)
pval_df <- melt(cor_res$P)
names(cor_df) <- c("Var1", "Var2", "Correlation")
names(pval_df) <- c("Var1", "Var2", "p_value")

# Merge correlation and p-value data
cor_full <- left_join(cor_df, pval_df, by = c("Var1", "Var2"))


# Keep lower triangle only (non-mirrored)
cor_full <- cor_full %>%
  filter(as.integer(factor(Var1, levels = unique(Var1))) > 
           as.integer(factor(Var2, levels = unique(Var2))))


ggplot(cor_full, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 4) +
  scale_fill_gradient2(low = "darkblue", high = "lightpink", mid = "white", midpoint = 0,
                       limit = c(-1, 1), name = "Pearson r") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Correlation Heatmap of Drought Response Ratios",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_full <- cor_full %>%
  mutate(sig = case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ ""
  ))

ggplot(cor_full, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(Correlation, 2), sig)), size = 4) +
  scale_fill_gradient2(low = "darkblue", high = "lightpink", mid = "white", midpoint = 0,
                       limit = c(-1, 1), name = "Pearson r") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Correlation Heatmap of Drought Response Ratios (Lower Triangle)",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add significance labels
cor_full <- cor_full %>%
  mutate(sig_label = case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ ""
  ),
  label_text = paste0(round(Correlation, 2), sig_label))

ggplot(cor_full, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label_text), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Correlation Heatmap of Drought Response Ratios",
       subtitle = "Significant correlations annotated (* p ≤ 0.05)",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







# Enhanced correlation matrix of Drought Response Ratios with Origins variables
library(dplyr)
library(tidyr)
library(Hmisc)
library(reshape2)
library(ggplot2)

# Define traits for response ratios
traits <- c("ShootWt", "Florets", "d15N", "d13C", "N", "P", "RLC_P")

# First, check what columns actually exist in Origins dataset
print("=== AVAILABLE COLUMNS IN ORIGINS DATASET ===")
print(colnames(Origins))

# Define potential continuous origins variables 
potential_origins_continuous <- c(
  "Origin Anthesis Date (DAP)*",
  "Anthesis date (days)",
  "Harvest date (days)", 
  "Total fresh weight (kg)",
  "Brix (maturity)",
  "Brix (milk)",
  "Dry weight (kg)", 
  "Stalk height (cm)",
  "Dry tons per acre",
  "ADF (% DM)",
  "NDF (% DM)",
  "NFC (% DM)",
  "Lignin (% DM)"
)

# Only keep origins variables that actually exist in the dataset
origins_continuous <- potential_origins_continuous[potential_origins_continuous %in% colnames(Origins)]
print("=== ORIGINS VARIABLES FOUND IN DATASET ===")
print(origins_continuous)

# Clean origins data and ensure Genotype format matches
Origins_clean <- Origins
Origins_clean$Genotype <- gsub("\\.0$", "", as.character(Origins_clean$Genotype))

# Convert continuous variables to numeric
for(col in origins_continuous) {
  if(col %in% colnames(Origins_clean)) {
    Origins_clean[[col]] <- as.numeric(as.character(Origins_clean[[col]]))
  }
}

# Combine macombo with origins data
macombo_with_origins <- macombo %>%
  left_join(Origins_clean, by = "Genotype")

# Compute drought response ratios and include origins data
rr_with_origins <- macombo_with_origins %>%
  group_by(Genotype) %>%
  summarise(
    # Calculate response ratios for traits
    across(all_of(traits), 
           ~ log(mean(.x[Treat == "Droughted"], na.rm = TRUE) / 
                   mean(.x[Treat == "Watered"], na.rm = TRUE)),
           .names = "rr_{.col}"),
    # Include origins variables (take first non-NA value per genotype)
    across(all_of(origins_continuous), ~ first(.x[!is.na(.x)]), .names = "orig_{.col}"),
    .groups = 'drop'
  )

# Get available origins variables with sufficient data
available_origins <- origins_continuous[
  sapply(origins_continuous, function(x) {
    col_name <- paste0("orig_", x)
    if(col_name %in% colnames(rr_with_origins)) {
      sum(!is.na(rr_with_origins[[col_name]])) > 5
    } else {
      FALSE
    }
  })
]

# Response ratio variable names
rr_vars <- paste0("rr_", traits)
orig_vars <- paste0("orig_", available_origins)

# Create final correlation dataset
cor_data_combined <- rr_with_origins %>%
  select(all_of(c(rr_vars, orig_vars))) %>%
  select_if(function(x) sum(is.finite(x)) > 5) %>%
  filter(if_all(everything(), is.finite))

print(paste("Number of genotypes with complete data:", nrow(cor_data_combined)))
print(paste("Response ratio variables:", length(rr_vars)))
print(paste("Origins variables included:", length(intersect(orig_vars, colnames(cor_data_combined)))))

# Function to create enhanced correlation heatmap
create_enhanced_heatmap <- function(data, title = "Correlation Matrix") {
  if(ncol(data) < 2) {
    print("Insufficient data for correlation matrix")
    return(NULL)
  }
  
  cor_result <- rcorr(as.matrix(data))
  
  cor_df <- melt(cor_result$r)
  pval_df <- melt(cor_result$P)
  names(cor_df) <- c("Var1", "Var2", "Correlation")
  names(pval_df) <- c("Var1", "Var2", "p_value")
  
  cor_full <- left_join(cor_df, pval_df, by = c("Var1", "Var2")) %>%
    mutate(
      sig_label = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**", 
        p_value <= 0.05 ~ "*",
        TRUE ~ ""
      ),
      label_text = ifelse(is.na(p_value), 
                          round(Correlation, 2),
                          paste0(round(Correlation, 2), sig_label)),
      # Categorize variables for better visualization
      var1_type = ifelse(grepl("^rr_", Var1), "Response Ratio", "Origin Trait"),
      var2_type = ifelse(grepl("^rr_", Var2), "Response Ratio", "Origin Trait")
    )
  
  # Clean variable names for display
  cor_full$Var1_clean <- gsub("^rr_|^orig_", "", cor_full$Var1)
  cor_full$Var2_clean <- gsub("^rr_|^orig_", "", cor_full$Var2)
  
  ggplot(cor_full, aes(x = Var2_clean, y = Var1_clean, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = label_text), size = 2.5) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
    theme_minimal() +
    coord_fixed() +
    labs(title = title,
         subtitle = "Significance: *** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05",
         x = "", y = "") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 12),
      legend.position = "right"
    )
}

# Create the enhanced correlation heatmap
if(nrow(cor_data_combined) > 5 && ncol(cor_data_combined) > 1) {
  enhanced_plot <- create_enhanced_heatmap(cor_data_combined, 
                                          "Drought Response Ratios + Origins Characteristics")
  print(enhanced_plot)
}

# Focus on cross-correlations between response ratios and origins only
if(length(intersect(rr_vars, colnames(cor_data_combined))) > 0 && 
   length(intersect(orig_vars, colnames(cor_data_combined))) > 0) {
  
  # Calculate correlation matrix
  full_cor_matrix <- cor(cor_data_combined, use = "complete.obs")
  
  # Extract response ratios vs origins correlations
  available_rr_vars <- intersect(rr_vars, colnames(cor_data_combined))
  available_orig_vars <- intersect(orig_vars, colnames(cor_data_combined))
  
  rr_origin_cors <- full_cor_matrix[available_rr_vars, available_orig_vars, drop = FALSE]
  
  # Melt for plotting
  rr_origin_df <- melt(rr_origin_cors)
  names(rr_origin_df) <- c("Response_Ratio", "Origin_Trait", "Correlation")
  
  # Clean names for display
  rr_origin_df$Response_Ratio_clean <- gsub("^rr_", "", rr_origin_df$Response_Ratio)
  rr_origin_df$Origin_Trait_clean <- gsub("^orig_", "", rr_origin_df$Origin_Trait)
  
  # Create focused heatmap
  focused_plot <- ggplot(rr_origin_df, aes(x = Origin_Trait_clean, y = Response_Ratio_clean, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
    theme_minimal() +
    coord_fixed() +
    labs(title = "Drought Response Ratios vs Origin Characteristics",
         x = "Origin Characteristics", y = "Drought Response Ratios") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(focused_plot)
  
  # Print strongest correlations
  print("\n=== STRONGEST RESPONSE RATIO - ORIGIN CORRELATIONS ===")
  strongest_cors <- rr_origin_df %>%
    arrange(desc(abs(Correlation))) %>%
    head(10)
  
  for(i in 1:min(5, nrow(strongest_cors))) {
    print(paste(strongest_cors$Response_Ratio_clean[i], "vs", 
                strongest_cors$Origin_Trait_clean[i], ":", 
                round(strongest_cors$Correlation[i], 3)))
  }
}

# Summary statistics
print(paste("\nSummary:"))
print(paste("- Genotypes with complete data:", nrow(cor_data_combined)))
print(paste("- Response ratio variables:", length(intersect(rr_vars, colnames(cor_data_combined)))))
print(paste("- Origins variables:", length(intersect(orig_vars, colnames(cor_data_combined)))))

if(exists("rr_origin_cors") && length(rr_origin_cors) > 0) {
  print(paste("- Mean absolute correlation (RR vs Origins):", 
              round(mean(abs(rr_origin_cors), na.rm = TRUE), 3)))
  print(paste("- Max absolute correlation:", 
              round(max(abs(rr_origin_cors), na.rm = TRUE), 3)))
}





#-------g2combo plots----

g2combo$TotPlantDryWt <- as.numeric(as.character(g2combo$TotPlantDryWt))

# Filter out NA and 0 values from amf
g2combosubset <- g2combo %>%
  filter(!is.na(TOT) & TOT != 0 & Treat != "Dead")

# Calculate the IQR and remove outliers
Q1 <- quantile(g2combosubset$RLC_Percent, 0.25)
Q3 <- quantile(g2combosubset$RLC_Percent, 0.75)
IQR <- Q3 - Q1

g2combosubset_no_outliers <- g2combosubset %>%
  filter(RLC_Percent >= (Q1 - 1.5 * IQR) & RLC_Percent <= (Q3 + 1.5 * IQR))

# Plot the filtered data
ggplot(data = g2combosubset_no_outliers, aes(x = Genotype, y = RLC_Percent)) +
  geom_boxplot(aes(fill = Genotype)) +  # Color boxplot by 'Genotype'
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate Genotype labels by 45 degrees
  ) + 
  scale_x_discrete(labels = function(x) gsub("\\.0$", "", x)) +  # Remove ".0" from labels
  stat_n_text() +  # Add sample size text inside the boxplot
  labs(y = "% Root Length Colonized", x = "Genotype",title = "G2C Root Length Colonized by Genotype")  # Y label and X label





#violin of Treat
g2combosubset_clean <- g2combosubset %>%
  filter(!is.na(Genotype), !is.na(Arb), !is.na(TOT)) %>%  # Remove rows with NA in 'Genotype', 'Arb', or 'TOT'
  mutate(Arb_perc = (Arb / TOT) * 100)  # Calculate percentage

# Plot
g2combosubset_clean %>%
  ggplot(aes(x = Genotype, y = Arb_perc, fill = Genotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with transparency
  stat_n_text() +  # Adds count labels on violin plots
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "% Root Length of Arbuscules", x = "")



g2combo %>%
  ggplot(aes(x = AM_Hyphae/TOT, y = Arb/TOT, fill = Treat, color = Treat)) +  # Color points and trendlines based on 'Treat'
  geom_jitter(width = 0.2, height = 0.2, size = 3, shape = 21, alpha = 0.7) +  # Customize jittered points
  geom_smooth(method = "lm", se = FALSE) +  # Add trendline with color based on 'Treat'
  theme_minimal() +  # Apply a minimal theme
  scale_fill_manual(values = c("Live" = "forestgreen", "Dead" = "darkred")) +  # Custom colors for points
  scale_color_manual(values = c("Live" = "forestgreen", "Dead" = "darkred")) +  # Custom colors for trendlines
  geom_text(
    aes(x = .1, y = .7, label = paste("Live: Slope =", round(slope_watered, 2), 
                                      "\nR² =", round(r2_watered, 2))),
    color = "lightblue", hjust = 0, size = 3
  ) +  # Annotate Watered slope and R²
  geom_text(
    aes(x = .1, y = 1, label = paste("Dead: Slope =", round(slope_droughted, 2), 
                                     "\nR² =", round(r2_droughted, 2))),
    color = "darkred", hjust = 0, size = 3
  )  # Annotate Droughted slope and R²






# Create the 'mgr' variable with the correct calculation
ds <- g2combo %>%
  group_by(Genotype, Treat) %>% # Grouping by Genotype and Treat
  mutate(mgr = log((TotPlantDryWt) / mean((TotPlantDryWt)[Treat == "Dead"], na.rm = TRUE))) %>% 
  ungroup() # Ungroup after mutate to avoid issues with further operations

# Remove outliers using IQR method
Q1 <- quantile(ds$mgr, 0.25, na.rm = TRUE)
Q3 <- quantile(ds$mgr, 0.75, na.rm = TRUE)
IQR_value <- IQR(ds$mgr, na.rm = TRUE)

ds_no_outliers <- ds %>%
  filter(mgr >= (Q1 - 1.5 * IQR_value) & mgr <= (Q3 + 1.5 * IQR_value))

# Arrange bars by mean of 'mgr' per 'Genotype'
ds_no_outliers <- ds_no_outliers %>%
  group_by(Genotype) %>%
  mutate(mean_mgr = mean(mgr, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = unique(Genotype[order(mean_mgr,decreasing= TRUE)])))

# Create the boxplot
ggplot(ds_no_outliers, aes(x = Genotype, y = mgr, fill = Genotype)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Horizontal line at y=0
  labs(title = "Total Weight Response Ratio", y = "Live Over Sterile Inoculation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), # Hide legend and remove x-axis label
        axis.text.x = element_text(angle = 45, hjust = 1)) # Angle x-axis labels by 45 degrees






#----------MAC 2022 Probe Data----

# Add a 'Source' column to differentiate between the two datasets
dsdrought$Source <- "Drought"
dswet$Source <- "Wet"


# Combine the datasets by stacking them together
combined_data <- bind_rows(dsdrought, dswet)

# Convert Timestamps to datetime (POSIXct) if not already done
# Assume Timestamps are currently in 'YYYY-MM-DD' format without times, add a dummy time if needed
dsdrought$Timestamps <- as.POSIXct(dsdrought$Timestamps, format = "%Y-%m-%d", tz = "UTC")
dswet$Timestamps <- as.POSIXct(dswet$Timestamps, format = "%Y-%m-%d", tz = "UTC")

# Reshape data from wide to long format for easy processing
combined_data_long <- combined_data %>%
  pivot_longer(cols = starts_with("Port"), 
               names_to = "Port", 
               values_to = "Value")

# Calculate the average and standard deviation for each Port, grouped by Source and Port
summary_data <- combined_data_long %>%
  group_by(Source, Port, Timestamps) %>%
  summarise(Average_Value = mean(Value, na.rm = TRUE),
            StdDev = sd(Value, na.rm = TRUE),
            .groups = "drop")

# Create the line graph with averaged values and error shadows (standard deviation)
ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Port, group = interaction(Source, Port))) +
  geom_line() +  # Create line plot
  geom_ribbon(aes(ymin = Average_Value - StdDev, ymax = Average_Value + StdDev, fill = Port), 
              alpha = 0.2) +  # Add error shadow (ribbon)
  facet_wrap(~ Source, scales = "free_y") +  # Facet by Port to have separate plots for each
  #grid.arrange( ncol = 1) +  # ncol = 1 arranges the plots in one column (top-bottom)
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +  # Remove X title
  labs(y = "m³/m³ Water Content at 25in")+  # Label y-axis and legend
  scale_x_datetime(
    breaks = scales::date_breaks("4 days"),  # Adjust for monthly breaks (change to "1 day", "2 weeks", etc.)
    labels = scales::date_format("%b %d"))  # Format with month-day-year hour:minute


# Create a new column to define whether a port is for Drought or Wet conditions
summary_data$Condition <- ifelse(as.numeric(gsub("Port", "", summary_data$Port)) <= 6, "Drought", "Wet")

# Define color scales for each condition
ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Port, group = interaction(Source, Port))) +
  geom_line(size = 1) +  # Make the lines thicker by setting the size (e.g., 1 or higher)
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +  # Stack the graphs vertically, one on top of the other
  theme_minimal() +  # Minimal theme for cleaner visualization
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none",  # No legend displayed
    panel.grid.major = element_line(size = 0.2, color = "gray80"),  # Lighter grid lines
    panel.grid.minor = element_blank()  # No minor grid lines
  ) +  
  labs(
    y = "m³/m³ Water Content at 25in",  # Y-axis label
    x = "Timestamp"  # X-axis label
  ) +
  scale_x_datetime(
    breaks = scales::date_breaks("2 days"),  # Adjust for breaks every 4 days
    labels = scales::date_format("%b %d")  # Format as month-day
  ) +
  # Define two separate gradients for color (one for Drought, one for Wet)
  scale_color_manual(
    values = c(
      # Red to Green for Drought (Port1 to Port6)
      setNames(colorRampPalette(c("darkred", "chartreuse3"))(6), paste0("Port", 1:6)),
      # Green to Blue for Wet (Port7 to Port12)
      setNames(colorRampPalette(c("darkorange", "darkblue"))(6), paste0("Port", 7:12))
    )
  )


ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Source, group = Source)) +
  geom_line() +  # Create line plot
  theme(axis.title.x = element_blank()) +  # Remove X title
  labs(y = "m³/m³ Water Content at 25in", color = "Treat") +  # Label y-axis and legend
  scale_x_datetime(
    breaks = scales::date_breaks("2 days"),  # Adjust for monthly breaks (change to "1 day", "2 weeks", etc.)
    labels = scales::date_format("%b %d")  # Format with month-day-year hour:minute
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )



#---------MAC 2022 EMH Graph-----

# Filter out 0 values and NA values
emh_filtered <- emh %>%
  filter(!is.na(AMFInDrySoil) & 
           AMFInDrySoil != 0 & 
           AMFInDrySoil != "#DIV/0!"&
         !is.na(DSEInDrySoil))

emh_filtered$AMFInDrySoil <- as.numeric(as.character(emh_filtered$AMFInDrySoil))

# Create the AMF boxplot
emh_filtered %>%
  ggplot(aes(x = Genotype, y = AMFInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "AMF in Dry Soil (m/g)")

# Create the DSE boxplot
emh_filtered %>%
  ggplot(aes(x = Genotype, y = DSEInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "DSE in Dry Soil (m/g)")

#violin of Treat
emh_filtered %>%
  ggplot(aes(x = Treat, y = AMFInDrySoil, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


#--------Mac 2022 Colonization Graphs------

# Filter out NA and 0 values from amf
amf_filtered <- amf %>%
  filter(!is.na(RLC_P) & RLC_P != 0) %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Ensure at least 5 non-NA values for each Genotype
  ungroup()

#Calculate the IQR and remove outliers
Q1 <- quantile(amf_filtered$RLC_P, 0.25)
Q3 <- quantile(amf_filtered$RLC_P, 0.75)
IQR <- Q3 - Q1

amf_filtered_no_outliers <- amf_filtered %>%
  filter(RLC_P >= (Q1 - 1.5 * IQR) & RLC_P <= (Q3 + 1.5 * IQR))

# Plot the filtered data
ggplot(data = amf_filtered_no_outliers,
       aes(x = Treat, y = RLC_P)) +
  geom_boxplot(aes(fill = Treat)) +  # color boxplot
  theme(axis.title.x = element_blank(),  # remove X title
        axis.text.x = element_blank()) +  # remove X axis labels
  facet_grid(~Genotype) +  # facet by Genotype
  stat_n_text() +  # Add sample size text
  labs(y = "% Root Length Colonized", x = "")

# box plot of Arbs
ggplot(data = amf_filtered,
       aes(x = Treat, y = (LSE / TOT)*100)) +
  geom_boxplot(aes(fill = Treat)) +  # color boxplot
  theme(axis.title.x = element_blank(),  # remove X title
        axis.text.x = element_blank()) +  # remove X axis labels
  facet_grid(~Genotype) +  # facet by Genotype
  stat_n_text() +  # Add sample size text
  labs(y = "% Root Length LSE", x = "")

#violin of Treat
amf_filtered %>%
  ggplot(aes(x = Treat, y = ((Arb+Ves_and_Arb)/TOT)*100, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "% Root Length of Arbuscules", x = "")


#------- MAC 2022 Colonization Power Analysis-----
library(pwr)
# Subset the two treatment groups
watered <- amf_filtered %>% filter(Treat == "Watered")
droughted <- amf_filtered %>% filter(Treat == "Droughted")

# Compute mean and SD for each
m1 <- mean(watered$RLC_P, na.rm = TRUE)
m2 <- mean(droughted$RLC_P, na.rm = TRUE)
sd1 <- sd(watered$RLC_P, na.rm = TRUE)
sd2 <- sd(droughted$RLC_P, na.rm = TRUE)

# Pooled SD
s_pooled <- sqrt((sd1^2 + sd2^2) / 2)

# Cohen's d
cohen_d <- abs(m1 - m2) / s_pooled
cohen_d
n1 <- sum(!is.na(watered$RLC_P))
n2 <- sum(!is.na(droughted$RLC_P))

# Power calculation
power_result <- pwr.t2n.test(n1 = n1,
                             n2 = n2,
                             d = cohen_d,
                             sig.level = 0.05,
                             alternative = "two.sided")
power_result

# Range of sample sizes
sample_sizes <- seq(5, 50, by = 1)

# Compute power curve
power_curve <- sapply(sample_sizes, function(n) {
  pwr.t2n.test(n1 = n, n2 = n, d = cohen_d, sig.level = 0.05)$power
})

# Plot
library(ggplot2)
power_df <- data.frame(SampleSize = sample_sizes, Power = power_curve)

ggplot(power_df, aes(x = SampleSize, y = Power)) +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(title = "Power Curve for Detecting Treatment Effect on RLC_P",
       x = "Sample Size per Group",
       y = "Statistical Power") +
  theme_minimal()


# Function to compute Cohen's d and power for each genotype
power_by_genotype <- amf_filtered %>%
  filter(Treat %in% c("Watered", "Droughted")) %>%
  group_by(Genotype) %>%
  summarise(
    m_watered = mean(RLC_P[Treat == "Watered"], na.rm = TRUE),
    m_droughted = mean(RLC_P[Treat == "Droughted"], na.rm = TRUE),
    sd_watered = sd(RLC_P[Treat == "Watered"], na.rm = TRUE),
    sd_droughted = sd(RLC_P[Treat == "Droughted"], na.rm = TRUE),
    n_watered = sum(!is.na(RLC_P[Treat == "Watered"])),
    n_droughted = sum(!is.na(RLC_P[Treat == "Droughted"])),
    .groups = "drop"
  ) %>%
  mutate(
    pooled_sd = sqrt((sd_watered^2 + sd_droughted^2) / 2),
    cohen_d = abs(m_watered - m_droughted) / pooled_sd,
    power = pmap_dbl(list(n1 = n_watered, n2 = n_droughted, d = cohen_d), ~ {
      if (is.na(..3) || ..1 < 2 || ..2 < 2) return(NA)
      pwr.t2n.test(n1 = ..1, n2 = ..2, d = ..3, sig.level = 0.05)$power
    })
  ) %>%
  select(Genotype, n_watered, n_droughted, cohen_d, power)

print(power_by_genotype)


# Required sample size for d = 0.8 (medium effect size)
pwr.t.test(d = 0.8, power = 0.8, sig.level = 0.05, type = "two.sample")


#------Ordination-----

# Define predictor variables
vars <- c("RLC_P", "ShootWt", "Florets", "d15N", "d13C", "N", "P")

# Filter complete cases for PCA
# Include Treat when subsetting
pca_data <- macombo %>%
  select(Genotype, Treat, all_of(vars)) %>%
  drop_na()

# PCA on trait matrix
trait_matrix <- pca_data %>% select(all_of(vars))
pca_result <- prcomp(trait_matrix, center = TRUE, scale. = TRUE)

# Add scores and metadata
scores <- as.data.frame(pca_result$x)
scores$Genotype <- pca_data$Genotype
scores$Treat <- pca_data$Treat


# Loadings
loadings <- as.data.frame(pca_result$rotation)
loadings$Trait <- rownames(loadings)

ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Genotype, shape = Treat), size = 3, alpha = 0.8) +
  
  # Trait vectors (loadings)
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               inherit.aes = FALSE,  # <- prevent ggplot from using global aes
               color = "black") +
  
  # Trait labels
  geom_text(data = loadings,
            aes(x = PC1 * 3.2, y = PC2 * 3.2, label = Trait),
            inherit.aes = FALSE,
            size = 4, color = "black") +
  
  labs(title = "PCA of Drought-Response Traits",
       x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "right")


#-----Ordination by genotype----
# Drought Response Significance Heatmap
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

# Define traits to analyze
traits <- c("ShootWt", "Florets", "d15N", "d13C", "N", "P", "RLC_P","")

# Function to perform t-test and calculate effect size for each genotype-trait combination
calculate_drought_effects <- function(data, traits) {
  
  results_list <- list()
  
  for(trait in traits) {
    if(!trait %in% colnames(data)) {
      next
    }
    
    trait_results <- data %>%
      filter(!is.na(.data[[trait]]) & !is.na(Treat) & !is.na(Genotype)) %>%
      filter(Treat %in% c("Watered", "Droughted")) %>%
      group_by(Genotype) %>%
      filter(n() >= 4) %>%  # Need at least 4 observations (2 per treatment)
      reframe(  # Changed from summarise to reframe
        trait = trait,
        n_watered = sum(Treat == "Watered"),
        n_droughted = sum(Treat == "Droughted"),
        mean_watered = mean(.data[[trait]][Treat == "Watered"], na.rm = TRUE),
        mean_droughted = mean(.data[[trait]][Treat == "Droughted"], na.rm = TRUE),
        sd_watered = sd(.data[[trait]][Treat == "Watered"], na.rm = TRUE),
        sd_droughted = sd(.data[[trait]][Treat == "Droughted"], na.rm = TRUE),
        # Effect size (Cohen's d)
        cohens_d = {
          pooled_sd <- sqrt(((n_watered - 1) * sd_watered^2 + (n_droughted - 1) * sd_droughted^2) / 
                              (n_watered + n_droughted - 2))
          if(pooled_sd > 0) {
            (mean_droughted - mean_watered) / pooled_sd
          } else {
            0
          }
        },
        # T-test - simplified approach
        p_value = {
          if(n_watered >= 2 & n_droughted >= 2) {
            tryCatch({
              watered_vals <- .data[[trait]][Treat == "Watered"]
              droughted_vals <- .data[[trait]][Treat == "Droughted"]
              test_result <- t.test(droughted_vals, watered_vals)
              test_result$p.value
            }, error = function(e) {
              1.0
            })
          } else {
            1.0
          }
        },
        t_stat = {
          if(n_watered >= 2 & n_droughted >= 2) {
            tryCatch({
              watered_vals <- .data[[trait]][Treat == "Watered"]
              droughted_vals <- .data[[trait]][Treat == "Droughted"]
              test_result <- t.test(droughted_vals, watered_vals)
              as.numeric(test_result$statistic)
            }, error = function(e) {
              0
            })
          } else {
            0
          }
        }
      ) %>%
      mutate(
        # Significance categories
        significance = case_when(
          p_value <= 0.001 ~ "***",
          p_value <= 0.01 ~ "**",
          p_value <= 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        # Effect direction and magnitude for coloring
        effect_score = case_when(
          p_value > 0.05 ~ 0,  # Non-significant = neutral
          cohens_d > 0 ~ pmin(abs(cohens_d) * (-log10(p_value)), 10),  # Positive = red (drought better)
          cohens_d < 0 ~ -pmin(abs(cohens_d) * (-log10(p_value)), 10), # Negative = blue (watered better)
          TRUE ~ 0
        ),
        # Response ratio for interpretation
        response_ratio = log(mean_droughted / mean_watered)
      )
    
    results_list[[trait]] <- trait_results
  }
  
  # Combine all results
  bind_rows(results_list)
}

# Calculate drought effects
drought_effects <- calculate_drought_effects(macombo, traits)

print("=== DROUGHT EFFECTS SUMMARY ===")
print(paste("Total genotype-trait combinations:", nrow(drought_effects)))
print(paste("Significant effects (p < 0.05):", sum(drought_effects$p_value < 0.05, na.rm = TRUE)))

# Create the heatmap
heatmap_plot <- drought_effects %>%
  ggplot(aes(x = trait, y = Genotype, fill = effect_score)) +
  geom_tile(color = "white", size = 0.5) +
  
  # Add significance stars
  geom_text(aes(label = significance), 
            color = "white", size = 3, fontface = "bold") +
  
  # Custom color scale: Blue for watered preference, Red for drought tolerance
  scale_fill_gradient2(
    low = "darkblue", 
    mid = "gray90", 
    high = "darkred",
    midpoint = 0,
    name = "Effect\nScore",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  labs(
    title = "Drought Response by Genotype and Trait",
    subtitle = "Red = Better under drought, Blue = Better when watered\nDarker colors = more significant differences",
    x = "Trait",
    y = "Genotype",
    caption = "*** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05, ns = not significant"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(heatmap_plot)

# Summary statistics
print("\n=== SUMMARY BY TRAIT ===")
summary_by_trait <- drought_effects %>%
  group_by(trait) %>%
  summarise(
    n_genotypes = n(),
    n_significant = sum(p_value < 0.05),
    pct_significant = round(100 * n_significant / n_genotypes, 1),
    n_drought_better = sum(effect_score > 0),
    n_watered_better = sum(effect_score < 0),
    mean_effect_size = round(mean(abs(cohens_d), na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(summary_by_trait)

print("\n=== SUMMARY BY GENOTYPE ===")
summary_by_genotype <- drought_effects %>%
  group_by(Genotype) %>%
  summarise(
    n_traits = n(),
    n_significant = sum(p_value < 0.05),
    n_drought_better = sum(effect_score > 0),
    n_watered_better = sum(effect_score < 0),
    overall_drought_tolerance = mean(effect_score, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(overall_drought_tolerance))

print(head(summary_by_genotype, 10))

# Create a simplified version focusing on significant effects only
significant_effects <- drought_effects %>%
  filter(p_value < 0.05)

if(nrow(significant_effects) > 0) {
  sig_heatmap <- significant_effects %>%
    ggplot(aes(x = trait, y = Genotype, fill = effect_score)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = significance), 
              color = "white", size = 3, fontface = "bold") +
    scale_fill_gradient2(
      low = "darkblue", 
      mid = "gray90", 
      high = "darkred",
      midpoint = 0,
      name = "Effect\nScore"
    ) +
    labs(
      title = "Significant Drought Response Effects Only",
      subtitle = "Red = Better under drought, Blue = Better when watered",
      x = "Trait",
      y = "Genotype"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  print(sig_heatmap)
}

# Export results for further analysis
print("\n=== TOP DROUGHT TOLERANT GENOTYPES ===")
top_drought_tolerant <- summary_by_genotype %>%
  filter(overall_drought_tolerance > 0) %>%
  head(5)
print(top_drought_tolerant)

print("\n=== TOP WATER-DEPENDENT GENOTYPES ===")
top_water_dependent <- summary_by_genotype %>%
  filter(overall_drought_tolerance < 0) %>%
  tail(5)
print(top_water_dependent)

#-------More PCA-----
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

# Define traits to analyze (including new ones)
traits <- c("ShootWt", "Florets", "d15N", "d13C", "N", "P", "RLC_P", "C", "CN_Ratio")

# Calculate CN_Ratio if C and N columns exist
if(exists("macombo")) {
  if("C" %in% colnames(macombo) && "N" %in% colnames(macombo)) {
    macombo$CN_Ratio <- macombo$C / macombo$N
    print("Created CN_Ratio using C/N")
  } else {
    print("Warning: C or N columns not found - CN_Ratio cannot be calculated")
  }
}

# Define trait-specific interpretations for drought tolerance
# TRUE = higher values under drought indicate better tolerance
# FALSE = lower values under drought indicate better tolerance
trait_interpretations <- list(
  "ShootWt" = TRUE,     # Higher shoot weight under drought = better tolerance
  "Florets" = TRUE,     # More florets under drought = better tolerance
  "d15N" = FALSE,       # Lower δ15N under drought = better N status = better tolerance
  "d13C" = TRUE,        # Higher δ13C under drought = better water use efficiency = better tolerance
  "N" = TRUE,           # Higher N under drought = better nutrient status = better tolerance
  "P" = TRUE,           # Higher P under drought = better nutrient status = better tolerance
  "RLC_P" = TRUE,       # Higher RLC_P under drought = better P utilization = better tolerance
  "C" = TRUE,           # Higher C under drought = better carbon accumulation = better tolerance
  "CN_Ratio" = FALSE    # Lower C/N under drought = better N status = better tolerance
)

# Function to perform t-test and calculate effect size for each genotype-trait combination
calculate_drought_effects <- function(data, traits) {
  
  results_list <- list()
  
  for(trait in traits) {
    if(!trait %in% colnames(data) || trait == "") {
      next
    }
    
    trait_results <- data %>%
      filter(!is.na(.data[[trait]]) & !is.na(Treat) & !is.na(Genotype)) %>%
      filter(Treat %in% c("Watered", "Droughted")) %>%
      group_by(Genotype) %>%
      filter(n() >= 4) %>%  # Need at least 4 observations (2 per treatment)
      reframe(  # Changed from summarise to reframe
        trait = trait,
        n_watered = sum(Treat == "Watered"),
        n_droughted = sum(Treat == "Droughted"),
        mean_watered = mean(.data[[trait]][Treat == "Watered"], na.rm = TRUE),
        mean_droughted = mean(.data[[trait]][Treat == "Droughted"], na.rm = TRUE),
        sd_watered = sd(.data[[trait]][Treat == "Watered"], na.rm = TRUE),
        sd_droughted = sd(.data[[trait]][Treat == "Droughted"], na.rm = TRUE),
        # Effect size (Cohen's d)
        cohens_d = {
          pooled_sd <- sqrt(((n_watered - 1) * sd_watered^2 + (n_droughted - 1) * sd_droughted^2) / 
                              (n_watered + n_droughted - 2))
          if(pooled_sd > 0) {
            (mean_droughted - mean_watered) / pooled_sd
          } else {
            0
          }
        },
        # T-test - simplified approach
        p_value = {
          if(n_watered >= 2 & n_droughted >= 2) {
            tryCatch({
              watered_vals <- .data[[trait]][Treat == "Watered"]
              droughted_vals <- .data[[trait]][Treat == "Droughted"]
              test_result <- t.test(droughted_vals, watered_vals)
              test_result$p.value
            }, error = function(e) {
              1.0
            })
          } else {
            1.0
          }
        },
        t_stat = {
          if(n_watered >= 2 & n_droughted >= 2) {
            tryCatch({
              watered_vals <- .data[[trait]][Treat == "Watered"]
              droughted_vals <- .data[[trait]][Treat == "Droughted"]
              test_result <- t.test(droughted_vals, watered_vals)
              as.numeric(test_result$statistic)
            }, error = function(e) {
              0
            })
          } else {
            0
          }
        }
      ) %>%
      mutate(
        # Significance categories
        significance = case_when(
          p_value <= 0.001 ~ "***",
          p_value <= 0.01 ~ "**",
          p_value <= 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        # CORRECTED: Trait-specific effect scoring for drought tolerance
        drought_tolerance_score = case_when(
          p_value > 0.05 ~ 0,  # Non-significant = neutral
          # For traits where higher values under drought = better tolerance
          trait %in% names(trait_interpretations)[sapply(trait_interpretations, isTRUE)] & 
            cohens_d > 0 ~ pmin(abs(cohens_d) * (-log10(p_value)), 10),  # Positive score
          # For traits where lower values under drought = better tolerance  
          trait %in% names(trait_interpretations)[sapply(trait_interpretations, isFALSE)] & 
            cohens_d < 0 ~ pmin(abs(cohens_d) * (-log10(p_value)), 10),  # Positive score
          # All other cases = negative score (poor tolerance)
          TRUE ~ -pmin(abs(cohens_d) * (-log10(p_value)), 10)
        ),
        # Keep original effect score for comparison
        original_effect_score = case_when(
          p_value > 0.05 ~ 0,
          cohens_d > 0 ~ pmin(abs(cohens_d) * (-log10(p_value)), 10),
          cohens_d < 0 ~ -pmin(abs(cohens_d) * (-log10(p_value)), 10),
          TRUE ~ 0
        ),
        # Response ratio for interpretation
        response_ratio = log(mean_droughted / mean_watered),
        # Interpretation of the response
        drought_response_interpretation = case_when(
          p_value > 0.05 ~ "No significant effect",
          trait %in% names(trait_interpretations)[sapply(trait_interpretations, isTRUE)] & 
            cohens_d > 0 ~ "Drought tolerant response",
          trait %in% names(trait_interpretations)[sapply(trait_interpretations, isFALSE)] & 
            cohens_d < 0 ~ "Drought tolerant response",
          TRUE ~ "Drought sensitive response"
        )
      )
    
    results_list[[trait]] <- trait_results
  }
  
  # Combine all results
  bind_rows(results_list)
}

# Calculate drought effects
drought_effects <- calculate_drought_effects(macombo, traits)

print("=== DROUGHT EFFECTS SUMMARY ===")
print(paste("Total genotype-trait combinations:", nrow(drought_effects)))
print(paste("Significant effects (p < 0.05):", sum(drought_effects$p_value < 0.05, na.rm = TRUE)))

# Calculate overall drought tolerance for ordering (using corrected scores)
drought_effects <- drought_effects %>%
  group_by(Genotype) %>%
  mutate(
    overall_drought_tolerance = mean(drought_tolerance_score, na.rm = TRUE),
    original_overall_tolerance = mean(original_effect_score, na.rm = TRUE)
  ) %>%
  ungroup()

# Create trait labels for better display
trait_labels <- c(
  "ShootWt" = "Shoot Weight",
  "Florets" = "Florets", 
  "d15N" = "δ15N",
  "d13C" = "δ13C",
  "N" = "N",
  "P" = "P", 
  "RLC_P" = "RLC_P",
  "C" = "C",
  "CN_Ratio" = "C/N"
)

# Print trait interpretations for reference
print("\n=== TRAIT INTERPRETATIONS FOR DROUGHT TOLERANCE ===")
for(trait in names(trait_interpretations)) {
  direction <- ifelse(trait_interpretations[[trait]], "Higher under drought = better", "Lower under drought = better")
  print(paste(trait, ":", direction))
}

# Create enhanced significant effects heatmap with CORRECTED drought tolerance scoring
significant_effects <- drought_effects %>%
  filter(p_value < 0.05)

if(nrow(significant_effects) > 0) {
  
  # Corrected heatmap using proper drought tolerance scores
  sig_heatmap_corrected <- significant_effects %>%
    mutate(trait_display = trait_labels[trait]) %>%
    ggplot(aes(x = factor(trait_display, levels = trait_labels[traits]), 
               y = reorder(Genotype, overall_drought_tolerance), fill = drought_tolerance_score)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = significance), 
              color = "white", size = 2.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#d73027",      # Red = drought sensitive
      mid = "#f7f7f7",      # Gray = neutral
      high = "#1a9850",     # Green = drought tolerant
      midpoint = 0,
      name = "Drought\nTolerance"
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = "Drought Tolerance by Genotype and Trait",
      subtitle = "Genotypes ordered by overall drought tolerance (most tolerant at top)",
      x = NULL,
      y = "Genotype",
      caption = "*** p≤0.001  ** p≤0.01  * p≤0.05  |  Green=drought tolerant, Red=drought sensitive"
    ) +
    theme_minimal() +
    theme(
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  print(sig_heatmap_corrected)
  
  # Comparison plot showing the difference
  comparison_plot <- significant_effects %>%
    select(Genotype, trait, original_effect_score, drought_tolerance_score) %>%
    pivot_longer(cols = c(original_effect_score, drought_tolerance_score), 
                 names_to = "score_type", values_to = "score") %>%
    mutate(
      score_type = case_when(
        score_type == "original_effect_score" ~ "Original (Incorrect)",
        score_type == "drought_tolerance_score" ~ "Corrected (Trait-Specific)"
      ),
      trait_display = trait_labels[trait]
    ) %>%
    ggplot(aes(x = factor(trait_display, levels = trait_labels[traits]), 
               y = Genotype, fill = score)) +
    geom_tile(color = "white", size = 0.2) +
    facet_wrap(~score_type, ncol = 2) +
    scale_fill_gradient2(
      low = "#d73027", mid = "#f7f7f7", high = "#1a9850",
      midpoint = 0, name = "Score"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title = "Comparison: Original vs Corrected Drought Tolerance Scoring",
      subtitle = "Notice how some genotypes change ranking with trait-specific interpretation",
      x = NULL, y = "Genotype",
      caption = "Green = drought tolerant response, Red = drought sensitive response"
    ) +
    theme_minimal() +
    theme(
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 8),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank()
    )
  
  print(comparison_plot)
  
} else {
  print("No significant effects found for heatmap")
}

# CORRECTED Summary statistics
print("\n=== CORRECTED SUMMARY BY TRAIT ===")
summary_by_trait <- drought_effects %>%
  group_by(trait) %>%
  summarise(
    n_genotypes = n(),
    n_significant = sum(p_value < 0.05),
    pct_significant = round(100 * n_significant / n_genotypes, 1),
    n_drought_tolerant = sum(drought_tolerance_score > 0),
    n_drought_sensitive = sum(drought_tolerance_score < 0),
    mean_effect_size = round(mean(abs(cohens_d), na.rm = TRUE), 2),
    trait_interpretation = ifelse(trait %in% names(trait_interpretations), 
                                  ifelse(trait_interpretations[[trait]], "Higher=better", "Lower=better"), 
                                  "Unknown"),
    .groups = 'drop'
  )

print(summary_by_trait)

print("\n=== CORRECTED SUMMARY BY GENOTYPE ===")
summary_by_genotype <- drought_effects %>%
  group_by(Genotype) %>%
  summarise(
    n_traits = n(),
    n_significant = sum(p_value < 0.05),
    n_drought_tolerant = sum(drought_tolerance_score > 0),
    n_drought_sensitive = sum(drought_tolerance_score < 0),
    overall_drought_tolerance = mean(drought_tolerance_score, na.rm = TRUE),
    original_tolerance_score = mean(original_effect_score, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(overall_drought_tolerance))

print(head(summary_by_genotype, 10))

print("\n=== TOP DROUGHT TOLERANT GENOTYPES (CORRECTED) ===")
top_drought_tolerant <- summary_by_genotype %>%
  filter(overall_drought_tolerance > 0) %>%
  head(5)
print(top_drought_tolerant)

print("\n=== TOP DROUGHT SENSITIVE GENOTYPES (CORRECTED) ===")
top_drought_sensitive <- summary_by_genotype %>%
  filter(overall_drought_tolerance < 0) %>%
  tail(5)
print(top_drought_sensitive)

print("\n=== RANKING COMPARISON ===")
ranking_comparison <- summary_by_genotype %>%
  select(Genotype, overall_drought_tolerance, original_tolerance_score) %>%
  mutate(
    corrected_rank = rank(-overall_drought_tolerance),
    original_rank = rank(-original_tolerance_score),
    rank_change = original_rank - corrected_rank
  ) %>%
  arrange(corrected_rank) %>%
  head(10)

print(ranking_comparison)
#
#
#
#
#-----PCA with Origins data NOT WORKING----
# Enhanced PCA including Origins data
library(dplyr)
library(ggplot2)
library(grid)

# First, check what columns actually exist in Origins dataset
print("=== AVAILABLE COLUMNS IN ORIGINS DATASET ===")
print(colnames(Origins))

# Define potential continuous origins variables 
potential_origins_continuous <- c(
  "Origin Anthesis Date (DAP)*",
  "Anthesis date (days)",
  "Harvest date (days)", 
  "Total fresh weight (kg)",
  "Brix (maturity)",
  "Brix (milk)",
  "Dry weight (kg)", 
  "Stalk height (cm)",
  "Dry tons per acre",
  "ADF (% DM)",
  "NDF (% DM)",
  "NFC (% DM)",
  "Lignin (% DM)"
)

# Only keep origins variables that actually exist in the dataset
origins_continuous <- potential_origins_continuous[potential_origins_continuous %in% colnames(Origins)]
print("=== ORIGINS VARIABLES FOUND IN DATASET ===")
print(origins_continuous)

# Clean the Origins data
Origins_clean <- Origins
Origins_clean$Genotype <- gsub("\\.0$", "", as.character(Origins_clean$Genotype))

# Convert continuous variables to numeric (only for existing columns)
for(col in origins_continuous) {
  Origins_clean[[col]] <- as.numeric(as.character(Origins_clean[[col]]))
}

# Define predictor variables (traits)
trait_vars <- c("RLC_P", "ShootWt", "Florets", "d15N", "d13C", "N", "P")

# Combine all variables for PCA
all_pca_vars <- c(trait_vars, origins_continuous)

# Check if Genotype column exists in both datasets
print("=== CHECKING GENOTYPE COLUMNS ===")
print(paste("Genotype in macombo:", "Genotype" %in% colnames(macombo)))
print(paste("Genotype in Origins:", "Genotype" %in% colnames(Origins_clean)))

if(!"Genotype" %in% colnames(macombo)) {
  print("Available columns in macombo:")
  print(colnames(macombo))
  stop("Genotype column not found in macombo dataset")
}

if(!"Genotype" %in% colnames(Origins_clean)) {
  print("Available columns in Origins_clean:")  
  print(colnames(Origins_clean))
  stop("Genotype column not found in Origins_clean dataset")
}

# Combine macombo with origins data
macombo_with_origins <- macombo %>%
  left_join(Origins_clean, by = "Genotype")

# Filter complete cases for PCA - calculate genotype means first to avoid treatment effects
pca_data_genotype_means <- macombo_with_origins %>%
  filter(!is.na(Genotype) & Genotype != "") %>%  # Remove missing genotypes
  group_by(Genotype) %>%
  summarise(
    # Calculate means for trait variables (only if they exist)
    across(all_of(trait_vars[trait_vars %in% colnames(.)]), ~ mean(.x, na.rm = TRUE)),
    # Take first non-NA value for origins (should be same for all reps of a genotype)
    across(all_of(origins_continuous[origins_continuous %in% colnames(.)]), ~ first(.x[!is.na(.x)])),
    .groups = 'drop'
  ) %>%
  # Keep only variables with sufficient data
  select_if(function(x) sum(is.finite(x)) > 5) %>%
  drop_na()

print("=== DATA PROCESSING RESULTS ===")
print(paste("Genotypes with complete data:", nrow(pca_data_genotype_means)))
print(paste("Columns in pca_data_genotype_means:", ncol(pca_data_genotype_means)))
if(nrow(pca_data_genotype_means) > 0) {
  print("Column names:")
  print(colnames(pca_data_genotype_means))
}
# Get final variable list (excluding Genotype)
if("Genotype" %in% colnames(pca_data_genotype_means)) {
  final_vars <- setdiff(colnames(pca_data_genotype_means), "Genotype")
} else {
  print("WARNING: Genotype column missing from processed data")
  final_vars <- colnames(pca_data_genotype_means)
}

print("=== FINAL VARIABLES IN PCA ===")
print(final_vars)
print(paste("Number of variables:", length(final_vars)))

# PCA on combined trait and origins matrix
if(length(final_vars) > 1 && nrow(pca_data_genotype_means) > 3) {
  
  pca_matrix <- pca_data_genotype_means %>% select(all_of(final_vars))
  pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
  
  # Add scores and metadata
  scores <- as.data.frame(pca_result$x)
  scores$Genotype <- pca_data_genotype_means$Genotype
  
  # Add treatment information back (for coloring/shaping)
  # Take mean position per genotype for visualization
  treatment_info <- macombo %>%
    filter(!is.na(Genotype) & !is.na(Treat) & Genotype != "" & Treat != "") %>%
    group_by(Genotype) %>%
    summarise(
      Primary_Treat = {
        treat_table <- table(Treat)
        if(length(treat_table) > 0) {
          names(sort(treat_table, decreasing = TRUE))[1]
        } else {
          "Unknown"
        }
      },
      .groups = 'drop'
    )
  
  scores <- scores %>%
    left_join(treatment_info, by = "Genotype") %>%
    mutate(Primary_Treat = ifelse(is.na(Primary_Treat), "Unknown", Primary_Treat))
  
  # Loadings with variable type identification
  loadings <- as.data.frame(pca_result$rotation)
  loadings$Variable <- rownames(loadings)
  loadings$Type <- ifelse(loadings$Variable %in% trait_vars, "Plant Trait", "Origin Trait")
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Create enhanced PCA plot
  pca_plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Genotype), size = 3, alpha = 0.8) +
    
    # Trait vectors (loadings) - color by type
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4, color = Type),
                 arrow = arrow(length = unit(0.15, "cm")),
                 inherit.aes = FALSE,
                 size = 0.8) +
    
    # Trait labels
    geom_text(data = loadings,
              aes(x = PC1 * 4.3, y = PC2 * 4.3, label = Variable, color = Type),
              inherit.aes = FALSE,
              size = 3, fontface = "bold") +
    
    labs(title = "PCA of Plant Traits + Origin Characteristics",
         subtitle = paste("Genotype means (n =", nrow(scores), "genotypes)"),
         x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")) +
    
    scale_color_manual(values = c("Plant Trait" = "darkblue", "Origin Trait" = "darkred"),
                       name = "Variable Type") +
    
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12))
  
  print(pca_plot)
  
  # Print PCA summary
  print("\n=== PCA SUMMARY ===")
  print(summary(pca_result))
  
  # Print loadings for first two PCs
  print("\n=== LOADINGS FOR PC1 AND PC2 ===")
  loadings_summary <- loadings %>%
    arrange(desc(abs(PC1))) %>%
    select(Variable, Type, PC1, PC2)
  print(loadings_summary)
  
  # Contribution plot - which variables contribute most to PC1 and PC2
  loadings_long <- loadings %>%
    select(Variable, Type, PC1, PC2) %>%
    pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "Loading")
  
  contribution_plot <- ggplot(loadings_long, aes(x = reorder(Variable, abs(Loading)), 
                                                 y = Loading, fill = Type)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~PC, scales = "free_x") +
    labs(title = "Variable Contributions to Principal Components",
         x = "Variables", y = "Loading") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(contribution_plot)
  
  # Alternative plot: separate colors for genotypes and shapes for treatment
  if("Primary_Treat" %in% colnames(scores) && sum(!is.na(scores$Primary_Treat)) > 0) {
    pca_plot_alt <- ggplot(scores, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = Genotype, shape = Primary_Treat), size = 3, alpha = 0.8) +
      
      # Trait vectors
      geom_segment(data = loadings,
                   aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4),
                   arrow = arrow(length = unit(0.15, "cm")),
                   inherit.aes = FALSE,
                   color = "black", size = 0.8) +
      
      # Trait labels with type-based colors
      geom_text(data = loadings,
                aes(x = PC1 * 4.3, y = PC2 * 4.3, label = Variable, color = Type),
                inherit.aes = FALSE,
                size = 3, fontface = "bold") +
      
      labs(title = "PCA of Plant Traits + Origin Characteristics",
           subtitle = "Genotype means with treatment information",
           x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")) +
      
      scale_color_manual(values = c("Plant Trait" = "darkblue", "Origin Trait" = "darkred"),
                         name = "Variable Type",
                         guide = guide_legend(override.aes = list(shape = NA, linetype = 1, size = 1))) +
      
      guides(
        shape = guide_legend("Treatment", override.aes = list(color = "black")),
        color = guide_legend("Genotype")
      ) +
      
      theme_minimal() +
      theme(legend.position = "right")
    
    print(pca_plot_alt)
  }
  
} else {
  print("Insufficient data for PCA analysis")
  print(paste("Variables available:", length(final_vars)))
  print(paste("Genotypes with complete data:", nrow(pca_data_genotype_means)))
}

#------t-Test------

t.test(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95, ...)


#-------Anova------

one.way <- aov(MHyphDrySoil ~ Genotype, data= emh_filtered)
summary(one.way)
plot(one.way)

two.way <- aov(MHyphDrySoil ~ Genotype * Treat, data= emh_filtered)
summary(two.way)
plot(two.way)

blocking <- aov(MHyphDrySoil ~ Genotype + Treat + Position, data= emh_filtered)
summary(blocking)
plot(blocking)


one.way <- aov(RLC_P ~ Genotype, data= amf_filtered)
summary(one.way)
plot(one.way)

two.way <- aov(RLC_P ~ Genotype * Treat, data= amf_filtered)
summary(two.way)
plot(two.way)

blocking <- aov(RLC_P ~ Genotype + Treat + Position, data= amf_filtered)
summary(blocking)
plot(blocking)


model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way", "interaction", "blocking")

aictab(model.set, modnames = model.names)

par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))


# 95% confidence interval

tukey.two.way<-TukeyHSD(two.way)

tukey.two.way

tukey.plot.aov<-aov(AP_AM_tot ~ TP:Genotype, data=MAC_P)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

Mean.AP.AM.tot.data <- MAC_P %>%
  group_by(TP, Genotype) %>%
  summarise(AM_tot = mean(AP_AM_tot))

Mean.AP.AM.tot.data$group <- c("a","b","b","c")

Mean.AP.AM.tot.data

two.way.plot <- ggplot(Mean.AP.AM.tot.data, 
                       aes(x = Genotype, y = AM_tot, group=TP)) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

two.way.plot <- two.way.plot +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=Mean.AP.AM.tot.data, aes(x=Genotype, y=AM_tot))

two.way.plot

#-----Legacy Plot----

dsLegacyfiltered <- dsLegacy %>%
  filter(!is.na(ShootWt))

# Calculate IQR and filter out outliers
Q1 <- quantile(dsLegacyfiltered$ShootWt, 0.25)
Q3 <- quantile(dsLegacyfiltered$ShootWt, 0.75)
IQR_value <- IQR(dsLegacyfiltered$ShootWt)

# Filter out outliers
dsLegacy_no_outliers <- dsLegacyfiltered %>%
  filter(ShootWt >= (Q1 - 1.5 * IQR_value) & ShootWt <= (Q3 + 1.5 * IQR_value))

# Create the boxplot without outliers
dsLegacy_no_outliers %>%
  ggplot(aes(x = GenoLetter, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Adds count labels on boxplots
  theme(axis.title.x = element_blank())  # Removes x-axis title

#violin of Treat
dsLegacy_no_outliers %>%
  ggplot(aes(x = Treat, y = ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

# Boxplot with Treat as fill
dsLegacy %>%
  ggplot(aes(x = GenoLetter, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank())  # Adds count labels on boxplots

#-------MAC Plots----

# Remove NAs
macfiltered <- mac %>%
  filter(!is.na(ShootWt))  %>%
  filter(!is.na(Genotype)) %>%
  filter(!is.na(Treat)) %>%
  filter(!is.na(Florets))

# Convert ShootWt to numeric if necessary
macfiltered$ShootWt <- as.numeric(as.character(macfiltered$ShootWt))
macfiltered$Florets <- as.numeric(as.character(macfiltered$Florets))

# Boxplot with Treat as fill
macfiltered %>%
  ggplot(aes(x = Genotype, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank(),  # Removes x-axis title
        axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels

macfiltered %>%
  ggplot(aes(x = Genotype, y = Florets, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank(),  # Removes x-axis title
        axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels


macfiltered %>%
  # Calculate the average Florets within each Genotype-Treat combination
  group_by(Genotype, Treat) %>%
  summarise(avg_Florets = mean(Florets, na.rm = TRUE)) %>%
  ungroup() %>%
  # Reorder Genotype based on the average Florets (median for ordering)
  mutate(Genotype = reorder(Genotype, avg_Florets, FUN = median)) %>%
  # Now join back to the original dataset to keep `Florets` for plotting
  left_join(macfiltered, by = c("Genotype", "Treat")) %>%
  # Create the boxplot with the reordered Genotype
  ggplot(aes(x = Genotype, y = Florets, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Adds count labels on boxplots
  theme(axis.title.x = element_blank(),  # Removes x-axis title
  axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels


# Violin plot 
macfiltered %>%
  ggplot(aes(x = Treat, y = ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

# Violin plot 
macfiltered %>%
  ggplot(aes(x = Treat, y = Florets, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

macfiltered %>%
  ggplot(aes(x = Treat, y = Florets, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


dmgrfiltered <- dmgr %>%
  filter(!is.na(ShootWt))  %>%
  filter(!is.na(Genotype)) %>%
  filter(!is.na(Treat)) %>%
  filter(!is.na(DryRootWt))

# Violin plot 
dmgrfiltered %>%
  ggplot(aes(x = Treat, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


dmgrfiltered %>%
  filter(NumOfPlants %in% c("1","2","3","4","5")) %>%
  ggplot(aes(x = GenoNum, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_boxplot()+
  stat_n_text()


dmgrfiltered %>%
  filter(NumOfPlants %in% c("1")) %>%
  ggplot(aes(x = GenoNum, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_boxplot()+
  stat_n_text()

dmgrfiltered$NumOfPlants <- as.numeric(as.character(dmgrfiltered$NumOfPlants))

# Compute median ShootWt for each NumOfPlants
order_df <- dmgrfiltered %>%
  filter(NumOfPlants %in% c("1", "2", "3")) %>%
  group_by(NumOfPlants) %>%
  summarize(median_totwt = median(DryRootWt+ShootWt, na.rm = TRUE)) %>%
  arrange(median_totwt) %>%
  mutate(NumOfPlants = factor(NumOfPlants, levels = NumOfPlants))

# Join the order_df with the original data to get the ordered levels
ds_ordered <- dmgrfiltered %>%
  filter(NumOfPlants %in% c("1", "2", "3")) %>%
  left_join(order_df, by = "NumOfPlants") %>%
  mutate(NumOfPlants = factor(NumOfPlants, levels = order_df$NumOfPlants))

# Plot
ggplot(dmgrfiltered, aes(x = GenoNum, y = DryRootWt+ShootWt/NumOfPlants, fill = NumOfPlants)) +
  geom_boxplot() +
  stat_summary(aes(label = ..count..), fun = "count", geom = "text", vjust = -0.5) +
  facet_grid(~NumOfPlants) +
  theme_minimal() +
  labs(title = "Boxplot of Shoot Weight by Genotype and Number of Plants",
       x = "Genotype Number",
       y = "Shoot Weight")


#------DMGR Plots------

ds <- dmgr %>%
  filter(NumOfPlants == "1")

# Create the 'mgr' variable with the correct calculation
ds <- dmgr %>%
  group_by(Genotype, Treat) %>% # Grouping by Genotype and Treat
  mutate(mgr = log((DryRootWt+ShootWt) / mean((DryRootWt+ShootWt)[Treat == "S"], na.rm = TRUE))) %>% 
  ungroup() # Ungroup after mutate to avoid issues with further operations

# Remove outliers using IQR method
Q1 <- quantile(ds$mgr, 0.25, na.rm = TRUE)
Q3 <- quantile(ds$mgr, 0.75, na.rm = TRUE)
IQR_value <- IQR(ds$mgr, na.rm = TRUE)

ds_no_outliers <- ds %>%
  filter(mgr >= (Q1 - 1.5 * IQR_value) & mgr <= (Q3 + 1.5 * IQR_value))

# Arrange bars by mean of 'mgr' per 'Genotype'
ds_no_outliers <- ds_no_outliers %>%
  group_by(Genotype) %>%
  mutate(mean_mgr = mean(mgr, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = unique(Genotype[order(mean_mgr,decreasing= TRUE)])))

# Create the boxplot
ggplot(ds_no_outliers, aes(x = Genotype, y = mgr, fill = Genotype)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Horizontal line at y=0
  labs(title = "Total Weight Response Ratio", y = "Live Over Sterile Inoculation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), # Hide legend and remove x-axis label
        axis.text.x = element_text(angle = 45, hjust = 1)) # Angle x-axis labels by 45 degrees


#-------- MAC 2022 amf ------

ggplot(data = g2combo,
       aes(x = Genotype, y =Olpidium/TOT))+
  geom_boxplot(aes(fill = Treat))+ # color boxplot
  theme(axis.title.x=element_blank()) #remove X title)
  #facet_grid(dmgr#facet_grid(~Genotype)+
  #stat_n_text()


amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup()

ggplot(data = amf_filtered,
  aes(x = Treat, y = RLC_P)) +
  geom_boxplot(aes(fill = Treat)) +  # color boxplot
  theme(axis.title.x = element_blank()) +  # remove X title
  facet_grid(~Genotype)+
  stat_n_text()



# Join the two datasets by a common column, e.g., 'Genotype'
amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup() %>%
  left_join(dmgr, by = "Genotype")  # Assuming 'Genotype' is the common column

ds <- amf_filtered %>%
  group_by(Genotype, Treat)%>% #We want to be doing MGR in groups by Genotype and Treat.
  mutate(mgr =  log(RLC_P / (mean(RLC_P[Treat == "Watered"], na.rm = TRUE))), 
         NA_real_)

ds %>%
  group_by(Genotype) %>%
  ggplot(aes(x = Genotype, 
             y = mgr, 
             fill = Genotype)) +
  geom_boxplot() +
  stat_n_text() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank())  # Hide legend and rotate x-axis text
#facet_grid(~ Treat, labeller = labeller(soil = soil_labels))

# Plot with faceting by NumOfPlants from the dmgr dataset
ggplot(data = ds,
       aes(x = Genotype, y = mgr)) +
  geom_boxplot(aes(fill = Genotype)) +  # Color boxplot
  theme(axis.title.x = element_blank()) +  # Remove X title
  facet_grid(~NumOfPlants)+ # Facet by NumOfPlants from the dmgr dataset
  stat_n_text() 


  
# Filter and prepare the data
amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 0) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup() %>%
  left_join(mac, by = "Genotype")  # Join with mac dataset on Genotype

# Calculate mgr (magnesium response ratio) by Genotype and Treat
ds <- amf_filtered %>%
  group_by(Genotype, Treat) %>%
  mutate(mgr = log(RLC_P / mean(RLC_P[Treat == "Watered"], na.rm = TRUE))) %>%
  ungroup()  # Remove grouping for further operations

# Boxplot for ShootWt with count labels (stat_n_text should work if it's from a custom package)
ggplot(ds, aes(x = Genotype, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Assuming stat_n_text() is part of a custom package or extension
  theme_bw() +  # Use the 'black and white' theme
  theme(legend.position = "none", axis.title.x = element_blank())  # Hide legend and remove x-axis title



