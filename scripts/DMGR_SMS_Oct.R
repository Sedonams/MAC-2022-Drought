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
library(scales)
library(viridis)
library(ggrepel)  # For smart label positioning
library(tibble)
library(lme4)
library(nlme)

#Set up working directory
getwd()

#Work path
setwd("H:/MAC-2022-Drought")
#laptop path
setwd("C:/Users/sedon/MAC-2022-Drought-1")

#ds <- read.csv("data/MAC22_cleaned.csv") #Using the cleaned up dataset. 
ds <- read.csv("data/dmgr_cleaned.csv") #Using the cleaned up dataset. 

View(ds)

theme_set(theme_bw())

#write.csv(ds, "H:/MAC-2022-Drought/dmgr_cleaned_again.csv", row.names = F)

# convert to date safely
date_cols <- c("seed_or_reseed_date", "sprout_date", "transplant_date", "harvest_date")
for (col in date_cols) {
  if (col %in% names(ds)) {
    # Specify the date format, e.g. "%Y-%m-%d" or "%m/%d/%Y" as needed
    ds[[col]] <- as.Date(as.character(ds[[col]]), format = "%Y-%m-%dT%H:%M:%OSZ")
  }
}

ds$sprout_or_transplant_to_harvest <- ds$sprout_or_transplant_to_harvest %>%
  as.integer()

str(ds)
colnames(ds)

dotchart(ds$total_biomass, group = factor(ds$treatment))
boxplot(dry_root_wt ~ sprout_or_transplant_to_harvest:treatment, data = ds)
boxplot(total_biomass_by_num ~ genotype:treatment, data = ds)
boxplot(total_biomass ~ genotype:treatment, data = ds)

# Linear model for total biomass by genotype and treatment
m <- lm(ds$total_biomass ~ ds$sprout_or_transplant_to_harvest)
#Type ‘m’ to look at the model coefficients.
m
str(m)
plot(m$residuals ~ m$fitted.values)

# another model for total biomass by number of plants 
n <- lm(ds$total_biomass ~ ds$genotype + ds$treatment)
n
str(n)
plot(n$residuals ~ n$fitted.values)



library(car)
ds$treatment <- as.factor(ds$treatment)
leveneTest(total_biomass ~ sprout_or_transplant_to_harvest, data = ds)
plot(m) # Example: plot residuals
plot(n) # Example: plot residuals
shapiro.test(m$residuals) 
shapiro.test(n$residuals)




# Calculate Cook's distance and create dataframe with observation numbers
cooks_df <- data.frame(
  obs = 1:length(cooks.distance(m)),
  distance = cooks.distance(m)
)

# Calculate threshold for influential points (commonly used is 4/n)
threshold <- 4/nrow(ds)

# Create Cook's distance plot
p_cooks <- ggplot(cooks_df, aes(x = obs, y = distance)) +
  geom_point(aes(color = distance > threshold), size = 3) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("steelblue", "red")) +
  labs(
    title = "Cook's Distance Plot",
    subtitle = "Red dashed line indicates threshold (4/n)",
    x = "Observation Number",
    y = "Cook's Distance"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Print plot
print(p_cooks)


# Identify influential points
influential <- which(cooks.distance(m) > threshold)
if(length(influential) > 0) {
  cat("\nInfluential observations (Cook's D >", round(threshold, 4), "):\n")
  print(ds[influential, c("genotype", "treatment", "rep", "total_biomass")])
}




############################
#trying somthing
MAC22_MLM_total_biomass <- lmer(total_biomass ~ treatment * genotype + sprout_or_transplant_to_harvest +(1|num_of_plants) + (1|rep), data = ds)
plot(residuals(MAC22_MLM_total_biomass)) # Example: plot residuals
plot(MAC22_MLM_total_biomass)
#


# Filter to groups with at least 2 observations, didnt work
ds_bartlett <- ds %>% 
  group_by(treatment, genotype) %>% 
  filter(n() >= 2) %>% 
  ungroup()
m <- bartlett.test(total_biomass ~ interaction(treatment, genotype), data = ds_bartlett)
# Bartlett test does not return residuals or fitted values, so plotting them is not applicable.
print(m)
####################################

hist (ds$dry_root_wt)
hist (ds$total_biomass_by_num)
hist (n$residuals)

hist (ds$total_biomass)
hist (m$residuals)


#Q-Q plot for m model 
qqnorm(m$residuals)
#It can be helpful to see a reference diagonal. Add the line by typing
qqline(m$residuals)

# Q-Q plot for n model 
qqnorm(n$residuals)
qqline(n$residuals)

# Shapiro-Wilk test for normality
shapiro.test(m$residuals)
shapiro.test(n$residuals)
# If p-value is less than 0.05, we reject the null hypothesis that the data is normally distributed. meaning data is not normal.

table(ds$total_biomass == 0)

plot(ds)
library(car)

# Check for multicollinearity using Variance Inflation Factor (VIF)
m <- lm(data=ds, total_biomass ~ treatment + genotype + num_of_plants + sprout_or_transplant_to_harvest )
vif(m)

plot(data=ds, sprout_or_transplant_to_harvest~total_biomass)
M.loess <- lm(sprout_or_transplant_to_harvest~total_biomass, data = ds)
Fit <- fitted(M.loess)
head(Fit)
lines(x = ds$total_biomass, y = Fit)

#plot(ds) not working
colnames(ds)

colnames(ds)

coplot(total_biomass ~ num_of_plants | treatment,
       data = ds,
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

#Check for autocorrelation in residuals if data are time-ordered
acf(m$residuals, main = "ACF of Model Residuals")

colnames(ds)

ggplot(ds, aes(x = shoot_wt, y = dry_root_wt)) +
  geom_point(size = 2, color = "deeppink") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw()

model1 <- lm(shoot_wt ~ dry_root_wt, data = ds)
par(mfrow = c(2, 2))
plot(model1)

#install.packages("performance")   # if not already installed
library(performance)
#install.packages("sjPlot")
library(sjPlot)
check_model(model1)

summary(model1)
tab_model(model1)
anova(model1)

null_model <- lm(shoot_wt ~ 1, data = ds)
anova(null_model, model1)

g = ggplot(ds, aes(y=shoot_wt, x=dry_root_wt)) + geom_point(size=2, color="deeppink") + theme_bw()
g= g+ geom_smooth(method="lm", color="black")
g

m=ggpredict(model = model1, terms = c("dry_root_wt")) 
plot(m, add.data=TRUE, color="deeppink")

ds$genotype <- as.factor(ds$genotype)
ggplot(ds, aes(x = treatment, y = total_biomass, fill = treatment)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

model2 <- lm(total_biomass ~ treatment, data = ds)
summary(model2)
tab_model(model2)

colnames(ds)


ds$shoot_weigher <- as.factor(ds$shoot_weigher)
ggplot(ds, aes(x = shoot_weigher, y = shoot_wt, fill = shoot_weigher)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

ds$root_weigher <- as.factor(ds$root_weigher)
ggplot(ds, aes(x = root_weigher, y = dry_root_wt, fill = root_weigher)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

model3 <- lm(shoot_wt ~ shoot_weigher, data = ds)
summary(model3)
# You should not trust these p-values to determine significance instead look at the F statistic for overall model fit
anova(model3)



model4 <- lm(total_biomass ~ sprout_or_transplant_to_harvest * treatment, data = ds)
summary(model4)
anova(model4)
# Type III ANOVA to assess main effects while accounting for interaction
Anova(model4, type = 3)
plot_model(model4, type = "pred", terms = c("sprout_or_transplant_to_harvest", "treatment"))

library(emmeans)
# Estimate slope for each sex group and test if slopes differ
slope_tests <- emtrends(model4, ~ treatment, var = "sprout_or_transplant_to_harvest")

# Compare slopes between sexes (test interaction)
pairs(slope_tests)
test(slope_tests)

model5 <- lm(total_biomass ~ treatment + genotype + num_of_plants + sprout_or_transplant_to_harvest, data = ds)
summary(model5)
Anova(model5, type=2)

#ggeffects not working
# Clean plot with confidence ribbon
ggplot(pred, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", size = 1.2) +
  geom_point(data = spar, aes(x = Wing, y = Weight), color = "deeppink", alpha = 0.4, size = 2) +
  labs(x = "Wing Length (mm)", y = "Predicted Weight (g)",
       title = "Model-Predicted Weight by Wing Length\n(Controlling for Tarsus and Head Length)") +
  theme_bw(base_size = 14)



ds1 <- ds %>% mutate(across(c(sprout_or_transplant_to_harvest, num_of_plants), list(s = scale)))
model7 <- lm(total_biomass ~ sprout_or_transplant_to_harvest_s + num_of_plants_s, data = ds1)
summary(model7)
car::vif(model7)









# Visualize pairwise relationships between key numeric variables



# Pairwise relationships between key numeric variables
numeric_vars <- ds %>%
  dplyr::select(dry_root_wt, shoot_wt, total_biomass, 
  sprout_or_transplant_to_harvest, num_of_plants) %>%
  # Remove any rows with all NA values
  filter(rowSums(!is.na(.)) > 0)

# Create enhanced pairs plot with modified settings
p <- ggpairs(numeric_vars,
        lower = list(continuous = wrap("points", alpha = 0.3)),
        diag = list(continuous = wrap("barDiag", bins = 30)),
        upper = list(continuous = wrap("cor", size = 3)),
        title = "Relationships Between Key Variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(size = 8)
  )

# Display plot
print(p)

# Save high-resolution version of the plot 
ggsave("C:/Users/sms689/Downloads/enhanced_pairs_plot.png", plot = p, width = 10, height = 10, dpi = 300)


######################################################

#----SHAP values for biomass predictors----

# Load required packages
library(xgboost)
library(shapr)
library(dplyr)

# Prepare data
  shap_data <- ds %>%
    select(
    total_biomass,  # response
    days_to_harvest,
    num_of_plants,
    dry_root_wt,
    shoot_wt
  ) %>%
  na.omit()

X <- shap_data[, setdiff(names(shap_data), "total_biomass")]
y <- shap_data$total_biomass
dtrain <- xgb.DMatrix(as.matrix(X), label = y)

# Train model
param <- list(objective = "reg:squarederror", max_depth = 3, eta = 0.1)
xgb_model <- xgb.train(param, dtrain, nrounds = 50)

# ---- SHAPR steps ----

# 1. Create explainer object
explainer <- setup(
  x_train = X,
  model = xgb_model,
  verbose = NULL
)

# 2. Get baseline prediction
phi0 <- mean(y)

# 3. Run explanation
shap_values <- explain(
  x_test = X,
  x_train = X,
  model = xgb_model,
  explainer = explainer,
  approach = "empirical",   
  phi0 = phi0   # <- key fix
)

# 4. Inspect
head(shap_values$shapley_values)

# Calculate SHAP values
shap_values <- explain(
  X,
  explainer,
  approach = "empirical",
  n_combinations = 1000  # Adjust based on computational resources
)

# Create SHAP summary plot
shap_long <- data.frame(
  feature = rep(colnames(X), each = nrow(X)),
  shap_value = unlist(shap_values$phi),
  feature_value = unlist(lapply(X, function(x) rep(x, length(colnames(X)))))
)

# Plot SHAP summary
ggplot(shap_long, aes(x = shap_value, y = reorder(feature, abs(shap_value)))) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_point(aes(color = feature_value), alpha = 0.5) +
  scale_color_viridis_c() +
  labs(
    title = "SHAP Values for Biomass Predictors",
    x = "SHAP Value (impact on model output)",
    y = "Feature",
    color = "Feature Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

# Save plot
ggsave("shap_summary_plot.png", width = 10, height = 6, dpi = 300)

# Print feature importance summary
importance_summary <- data.frame(
  feature = colnames(X),
  mean_abs_shap = colMeans(abs(shap_values$phi))
) %>%
  arrange(desc(mean_abs_shap))

print(importance_summary)




####################################################





















#------Making RLC's for the colonization data----

ds <- ds %>%
  mutate(arb_ves_and_arb = rowSums(select(., arb, ves_and_arb), na.rm = TRUE)) %>%
  mutate(vesicle_or_spore_ves_and_arb = rowSums(select(., vesicle_or_spore, ves_and_arb), na.rm = TRUE)) %>%
  mutate(am_hyphae_dse_and_am = rowSums(select(., am_hyphae, dse_and_am), na.rm = TRUE)) %>%
  mutate(dse_dse_and_am = rowSums(select(., dse, dse_and_am), na.rm = TRUE))

# Count non-zero values in each column to confirm successful combination
sum(ds$arb_ves_and_arb != 0, na.rm = TRUE)
sum(ds$arb != 0, na.rm = TRUE)

columns_to_convert <- c("am_hyphae_dse_and_am","arb_ves_and_arb", "dse_dse_and_am", "vesicle_or_spore_ves_and_arb", "lse", "coil","fine_endo")

ds <- ds %>%
  mutate(across(all_of(columns_to_convert), 
                ~ if_else(!is.na(.x) & !is.na(tot), 
                          (.x / tot) * 100, 
                          NA_real_),
                .names = "{.col}_rlc"))

ds <- ds %>%
  rename(amf_rlc = rlc_p)

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
    axis.text.x.top = element_text(angle = 15, hjust = .2, vjust = 0, size = 10, face = "bold"),
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

p

# Optional: Save high-quality version
ggsave("drought_response_heatmap_simple.png", plot = p, width = 14, height = 8, dpi = 300, bg = "white")



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
    caption = "Significance: *** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05"
  ) +
  theme_minimal() +
  theme(
    # Text styling
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
    
    # Axis styling
    axis.text.x.top = element_text(angle = 15, hjust = .2, vjust = 0, size = 10, face = "bold"),
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

p

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
    title = "Principal Component Analysis Treatment X Genotype (Preliminary AMF Colonization Genotypes)",
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
    panel.grid.minor.x = element_blank(),
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
  filter(rowSums(is.na(select(., all_of(all_vars)))) <= length(all_vars) * 0.65) %>%
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
    title = "Principal Component Analysis Treatment X Genotype (Preliminary AMF Colonization Genotypes)",
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
    panel.grid.minor.x = element_blank(),
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

#---Boxplots of colonization----
amf_filtered <- ds %>%
  filter(!is.na(amf_rlc) & amf_rlc != 0) %>%
  group_by(genotype) %>%
  filter(sum(!is.na(amf_rlc)) >= 5) %>%  # Ensure at least 5 non-NA values for each Genotype
  ungroup()

#Calculate the IQR and remove outliers
Q1 <- quantile(amf_filtered$amf_rlc, 0.25)
Q3 <- quantile(amf_filtered$amf_rlc, 0.75)
IQR <- Q3 - Q1

amf_filtered_no_outliers <- amf_filtered %>%
  filter(amf_rlc >= (Q1 - 1.5 * IQR) & amf_rlc <= (Q3 + 1.5 * IQR))

# Plot the filtered data
p <- ggplot(amf_filtered_no_outliers, aes(x = treatment, y = amf_rlc)) +
  geom_boxplot(aes(fill = treatment), 
               alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = treatment), 
              width = 0.2, size = 1.8, alpha = 1) +
  facet_grid(~genotype, scales = "free_x", space = "free_x") +  # keep treatments aligned
  stat_n_text() +
  labs(y = "AMF Root Length Colonized (%)", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),   # spacing between panels
    panel.border = element_rect(color = "black", fill = NA)  # border around each facet
  )
p


#------Correlating Colonization to Drought Tolerance per Genotype----

# Define colonization variables and drought response variables
colonization_vars <- c("amf_rlc", "am_hyphae_dse_and_am_rlc", "arb_ves_and_arb_rlc", 
                       "dse_dse_and_am_rlc", "vesicle_or_spore_ves_and_arb_rlc", 
                       "lse_rlc", "coil_rlc", "olpidium_rlc", "mold_rlc", 
                       "plasmodiophorids_rlc", "dot_line_rlc", "non_am_rlc", 
                       "fine_endo_rlc", "amf_in_dry_soil", "dse_in_dry_soil")

drought_response_vars <- c("shoot_wt", "florets", "d13c", "c", "d15n", "n", "p", 
                           "cn_ratio", "np_ratio")

# Step 1: Calculate drought effect (standardized) for each genotype and variable
drought_effects <- ds %>%
  select(genotype, treatment, all_of(c(colonization_vars, drought_response_vars))) %>%
  # Scale within genotype
  group_by(genotype) %>%
  mutate(across(all_of(c(colonization_vars, drought_response_vars)), 
                ~ scale(.)[,1], .names = "{.col}_scaled")) %>%
  ungroup() %>%
  # Calculate drought effect (Droughted - Watered) for each variable
  select(genotype, treatment, ends_with("_scaled")) %>%
  pivot_longer(cols = ends_with("_scaled"), names_to = "variable", values_to = "scaled_value") %>%
  filter(!is.na(scaled_value)) %>%
  # Remove the "_scaled" suffix
  mutate(variable = str_remove(variable, "_scaled")) %>%
  # Calculate mean by treatment
  group_by(genotype, variable, treatment) %>%
  summarise(mean_scaled = mean(scaled_value, na.rm = TRUE), .groups = "drop") %>%
  # Pivot to get drought effect
  pivot_wider(names_from = treatment, values_from = mean_scaled) %>%
  mutate(drought_effect = Droughted - Watered) %>%
  filter(!is.na(drought_effect)) %>%
  select(genotype, variable, drought_effect)

# Step 2: Separate colonization and drought response effects
colonization_effects <- drought_effects %>%
  filter(variable %in% colonization_vars) %>%
  pivot_wider(names_from = variable, values_from = drought_effect, names_prefix = "col_")

drought_response_effects <- drought_effects %>%
  filter(variable %in% drought_response_vars) %>%
  pivot_wider(names_from = variable, values_from = drought_effect, names_prefix = "resp_")

# Step 3: Merge the datasets
correlation_data <- inner_join(colonization_effects, drought_response_effects, by = "genotype")

# Step 4: Calculate correlations for each genotype (if enough data)
# First, let's see which genotypes have enough data
genotype_data_summary <- correlation_data %>%
  group_by(genotype) %>%
  summarise(
    n_col_vars = sum(!is.na(select(., starts_with("col_")))),
    n_resp_vars = sum(!is.na(select(., starts_with("resp_")))),
    .groups = "drop"
  ) %>%
  filter(n_col_vars >= 3 & n_resp_vars >= 3)  # Need at least 3 variables of each type

print("Genotypes with sufficient data for correlation analysis:")
print(genotype_data_summary)

# Step 5: Calculate correlations across all genotypes (pooled analysis)
col_data <- correlation_data %>% select(starts_with("col_"))
resp_data <- correlation_data %>% select(starts_with("resp_"))

# Remove columns with all NAs
col_data_clean <- col_data[, colSums(!is.na(col_data)) > 0]
resp_data_clean <- resp_data[, colSums(!is.na(resp_data)) > 0]

# Calculate correlation matrix
if(ncol(col_data_clean) > 0 & ncol(resp_data_clean) > 0) {
  correlation_matrix <- cor(col_data_clean, resp_data_clean, use = "pairwise.complete.obs")
  
  # Convert to long format for plotting
  cor_long <- correlation_matrix %>%
    as.data.frame() %>%
    rownames_to_column("colonization_var") %>%
    pivot_longer(-colonization_var, names_to = "response_var", values_to = "correlation") %>%
    filter(!is.na(correlation))
  
  # Clean up variable names
  cor_long <- cor_long %>%
    mutate(
      colonization_var = str_remove(colonization_var, "col_"),
      response_var = str_remove(response_var, "resp_")
    )
  
  # Create prettier names
  pretty_col_names <- c(
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
    "dse_in_dry_soil" = "DSE in Dry Soil"
  )
  
  pretty_resp_names <- c(
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
  
  cor_long <- cor_long %>%
    mutate(
      colonization_pretty = ifelse(colonization_var %in% names(pretty_col_names),
                                   pretty_col_names[colonization_var],
                                   colonization_var),
      response_pretty = ifelse(response_var %in% names(pretty_resp_names),
                               pretty_resp_names[response_var],
                               response_var)
    )
  
  # Create correlation heatmap
  p_correlation <- ggplot(cor_long, aes(x = response_pretty, y = colonization_pretty, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", correlation)), 
              color = ifelse(abs(cor_long$correlation) > 0.5, "white", "black"), 
              size = 3) +
    scale_fill_gradient2(
      low = "#313695", mid = "white", high = "#A50026",
      midpoint = 0,
      name = "Correlation\n(r)",
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5)
    ) +
    labs(
      title = "Correlation: Colonization Drought Response vs Plant Drought Response",
      subtitle = "How fungal colonization changes correlate with plant trait changes under drought",
      x = "Plant Response Variables (Drought Effect)",
      y = "Colonization Variables (Drought Effect)",
      caption = "Values show Pearson correlation coefficients between drought effects"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.caption = element_text(size = 9, color = "gray50")
    )
  
  p_correlation
  ggsave("colonization_drought_correlation_heatmap.png", plot = p_correlation, 
         width = 12, height = 10, dpi = 300, bg = "white")
}

# Step 6: Individual genotype analysis (for genotypes with sufficient data)
individual_correlations <- list()

for(geno in genotype_data_summary$genotype) {
  geno_data <- correlation_data %>% filter(genotype == geno)
  
  geno_col <- geno_data %>% select(starts_with("col_"))
  geno_resp <- geno_data %>% select(starts_with("resp_"))
  
  # Remove columns with all NAs for this genotype
  geno_col_clean <- geno_col[, colSums(!is.na(geno_col)) > 0]
  geno_resp_clean <- geno_resp[, colSums(!is.na(geno_resp)) > 0]
  
  if(ncol(geno_col_clean) >= 2 & ncol(geno_resp_clean) >= 2) {
    tryCatch({
      geno_cor <- cor(geno_col_clean, geno_resp_clean, use = "pairwise.complete.obs")
      individual_correlations[[geno]] <- geno_cor
    }, error = function(e) {
      cat("Error calculating correlation for genotype", geno, ":", e$message, "\n")
    })
  }
}

# Step 7: Summary statistics
cat("\n=== CORRELATION ANALYSIS SUMMARY ===\n")
cat("Number of genotypes analyzed:", nrow(correlation_data), "\n")
cat("Colonization variables available:", paste(names(col_data_clean), collapse = ", "), "\n")
cat("Response variables available:", paste(names(resp_data_clean), collapse = ", "), "\n")

if(exists("correlation_matrix")) {
  # Find strongest correlations
  strong_correlations <- cor_long %>%
    filter(abs(correlation) > 0.3) %>%
    arrange(desc(abs(correlation)))
  
  cat("\nStrongest correlations (|r| > 0.3):\n")
  if(nrow(strong_correlations) > 0) {
    for(i in 1:min(10, nrow(strong_correlations))) {
      cat(sprintf("%s <-> %s: r = %.3f\n", 
                  strong_correlations$colonization_pretty[i],
                  strong_correlations$response_pretty[i],
                  strong_correlations$correlation[i]))
    }
  } else {
    cat("No correlations with |r| > 0.3 found\n")
  }
}

# Step 8: Create a scatter plot for the strongest correlation
if(exists("strong_correlations") && nrow(strong_correlations) > 0) {
  # Get the strongest correlation
  strongest <- strong_correlations[1, ]
  
  # Get the actual data for plotting
  col_var <- paste0("col_", strongest$colonization_var)
  resp_var <- paste0("resp_", strongest$response_var)
  
  scatter_data <- correlation_data %>%
    select(genotype, all_of(c(col_var, resp_var))) %>%
    rename(x_var = all_of(col_var), y_var = all_of(resp_var)) %>%
    filter(!is.na(x_var) & !is.na(y_var))
  
  if(nrow(scatter_data) > 3) {
    p_scatter <- ggplot(scatter_data, aes(x = x_var, y = y_var)) +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
      geom_text_repel(aes(label = genotype), size = 3, alpha = 0.8) +
      labs(
        title = "Strongest Correlation: Colonization vs Drought Response",
        subtitle = paste0(strongest$colonization_pretty, " vs ", strongest$response_pretty, 
                          " (r = ", sprintf("%.3f", strongest$correlation), ")"),
        x = paste(strongest$colonization_pretty, "Drought Effect (standardized)"),
        y = paste(strongest$response_pretty, "Drought Effect (standardized)"),
        caption = "Each point represents one genotype's drought response"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.title = element_text(size = 11, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      )
    
    p_scatter
    ggsave("strongest_correlation_scatter.png", plot = p_scatter, 
           width = 10, height = 8, dpi = 300, bg = "white")
  }
}

# Step 9: Export correlation results
if(exists("cor_long")) {
  write.csv(cor_long, "colonization_drought_correlations.csv", row.names = FALSE)
  cat("\nCorrelation results exported to 'colonization_drought_correlations.csv'\n")
}







#---Preliminary Linear Model
str(ds)

# Look at the distribution of samples for each factor
table(ds[ , c("treatment", "genotype")])

#remove na's for the variables of interest
ds <- ds %>%
  filter(!is.na(sprout_or_transplant_to_harvest
) & !is.na(total_biomass))

# Look at the distribution of continuous variables:
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
hist(ds$sprout_or_transplant_to_harvest, xlab = "sprout_or_transplant_to_harvest", 
main = "")
hist(ds$total_biomass, xlab = "total_biomass", main = "")

# Standardized total_biomass, "by hand"
ds$Z_tot <- (ds$total_biomass - mean(ds$total_biomass)) / 
                      sd(ds$total_biomass)

hist(ds$Z_tot, xlab = "Z_tot", main = "")


# Standardized age, with the function scale
ds$Z_sprout_or_transplant_to_harvest <- scale(ds$sprout_or_transplant_to_harvest)

hist(ds$Z_sprout_or_transplant_to_harvest, xlab = "Z_sprout_or_transplant_to_harvest", main = "")

## Create a linear model without random effects
lm.test <- lm(Z_tot ~ Z_sprout_or_transplant_to_harvest, data = ds)
summary(lm.test)

## plot the data with ggplot
prelim_plot <- ggplot(ds, aes(x = Z_sprout_or_transplant_to_harvest, y = Z_tot, color = treatment)) +
  geom_point() + theme_bw()+
  geom_smooth(method = "lm", se=FALSE)+
    facet_wrap(~genotype)
prelim_plot

prelim_plot <- ggplot(ds, aes(x = sprout_or_transplant_to_harvest, y = total_biomass, color = treatment)) +
  geom_point() + theme_bw()+
  geom_smooth(method = "lm", se=FALSE)
prelim_plot

prelim_plot <- ggplot(ds, aes(x = num_of_plants, y = total_biomass, color = treatment)) +
  geom_point() + theme_bw()+
  geom_smooth(method = "lm", se=FALSE)+ 
  facet_wrap(~genotype)
prelim_plot


## Calculate residuals of this linear model
lm.test.resid <- rstandard(lm.test)

M.gls <- gls(Z_tot ~ Z_sprout_or_transplant_to_harvest, data = ds)

plot(lm.test.resid ~ as.factor(ds$genotype),
     xlab = "genotype", ylab = "Standardized residuals")+ abline(0, 0, lty = 2)


library(nlme)
M1.lme= lme(Z_tot ~ Z_sprout_or_transplant_to_harvest, random=~1|treatment, 
     data = ds, method="REML")

# Now add random effects
M1.lmer <- lmer(Z_tot ~ Z_sprout_or_transplant_to_harvest + (1 | genotype), data = ds)
              
summary(M1.lmer)


anova(M.gls, M1.lme)


# Plot predicted values vs residual values
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(M1.lme) ~ fitted(M1.lme), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

# In order to check the independence of the model residuals
# we need to plot residuals vs each covariate of the model

plot(resid(M1.lme) ~ ds$Z_tot, xlab = "Z_tot", ylab = "Normalized residuals")+abline(h = 0, lty = 2)

boxplot(resid(M1.lme) ~ treatment, data = ds, xlab = "treatment", ylab = "Normalized residuals") 

boxplot(resid(M1.lme) ~ genotype, data = ds, xlab = "genotype",
    ylab = "Normalized residuals")

hist(resid(M1.lme))

re <- ranef(M1.lmer, condVar = TRUE)
re_df <- as.data.frame(re$treatment)
ggplot(re_df, aes(sample = `(Intercept)`)) +
  stat_qq() + stat_qq_line() +
  labs(title = "QQ-plot of treatment random intercepts")

re <- ranef(M1.lmer, condVar = TRUE)
re_df <- as.data.frame(re$genotype)
ggplot(re_df, aes(sample = `(Intercept)`)) +
  stat_qq() + stat_qq_line() +
  labs(title = "QQ-plot of treatment random intercepts")

library(performance)
check_model(M1.lme)
check_model(M1.lmer)


isSingular(M1.lmer)

summary(M1.lme)


M1.lme= lme(Z_sprout_or_transplant_to_harvest ~ Z_tot, random=~1|treatment, 
     data = ds, method="ML")

M2.lme= lme(Z_sprout_or_transplant_to_harvest ~ 1, random=~1|treatment, 
     data = ds, method="ML")

anova(M1.lme, M2.lme)


M1.lme= lme(Z_sprout_or_transplant_to_harvest ~ Z_tot, random=~1|genotype, 
     data = ds, method="REML")
summary(M1.lme)

# ggplot theme

fig <- theme_bw() + theme(panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), panel.background = element_blank()) +
    theme(strip.background = element_blank(), strip.text.y = element_text()) +
    theme(legend.background = element_blank()) + theme(legend.key = element_blank()) +
    theme(panel.border = element_rect(colour = "black", fill = NA))

library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
pred.mm <- ggpredict(M1.lme, terms = c("Z_tot"))  # this gives overall predictions for the model
head(pred.mm)

# Plot the predictions 

ggplot(pred.mm) + 
   geom_line(aes(x = x, y = predicted)) +          # slope
   geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
               fill = "lightgrey", alpha = 0.5) +  # error band
   geom_point(data = ds,                      # adding the raw data (scaled values)
              aes(x = Z_tot, y = Z_sprout_or_transplant_to_harvest, colour = treatment)) + 
   labs(x = "tot wt", y = "sprout_or_transplant_to_harvest", 
        title = "All data") + 
   fig

install.packages("MuMIn")
library(MuMIn)
r.squaredGLMM(M1.lme)


M.randslope <- lme(Z_sprout_or_transplant_to_harvest ~ Z_tot,
                   random = ~ Z_tot | genotype,
                   data = ds,
                   method = "REML")

anova(M1.lme, M.randslope)


m3 <- lmer((total_biomass) ~ (treatment) + (genotype) + (1 | sprout_or_transplant_to_harvest) + (1| rep) , data = ds)
summary(m3)


Mbea <- lmer((total_biomass) ~ (treatment) + (genotype) +  (1| rep) + (1 | sprout_or_transplant_to_harvest) + (1 | purpleness), data = ds)
summary(Mbea) 

names(ds)

p <- ggplot(ds, aes(x = total_biomass, y = numeric_pt)) +
  geom_boxplot(aes(fill = genotype), 
               alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = genotype), 
              width = 0.2, size = 1.8, alpha = 1) +
  facet_grid(~treatment, scales = "free_x", space = "free_x") +  # keep treatments aligned
  stat_n_text() +
  #labs(y = "AMF Root Length Colonized (%)", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),   # spacing between panels
    panel.border = element_rect(color = "black", fill = NA)  # border around each facet
  )
p

p <- ggplot(ds, aes(x = purpleness, y = numeric_pt)) +
  geom_point(aes(color = genotype), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  facet_wrap(~treatment) +
  #labs(x = "Total Biomass", y = "Days to Harvest") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
p

view(ds)


# Option A: treat num_of_plants and sprout_or_transplant_to_harvest as random effects
m1 <- lmer(total_biomass ~ treatment * genotype +
             (1 | rep) +
             (1 | num_of_plants) +
             (1 | sprout_or_transplant_to_harvest),
           data = ds)

# Option B: treat them as covariates (fixed effects)
m2 <- lmer(dry_root_wt ~ treatment * genotype
             + sprout_or_transplant_to_harvest +
             (1 | rep),
           data = ds)

# Check models
summary(m1)
summary(m2)

# Compare fit
anova(m1, m2)

library(performance)

check_model(m2)

library(ggplot2)

ggplot(ds, aes(x = sprout_or_transplant_to_harvest, y = dry_root_wt)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    x = "Days from sprout/transplant to harvest",
    y = "Dry root weight (g)",
    title = "Dry root weight vs. harvest age"
  ) +
  theme_minimal()
  
  ggplot(ds, aes(x = sprout_or_transplant_to_harvest, y = dry_root_wt, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

dredge(  m2 )

na.omit(ds)
library(MuMIn)
library(lme4)
library(lmerTest)
library(nlme)

GGally::ggpairs(ds[, c("total_biomass", "dry_root_wt", "sprout_or_transplant_to_harvest", "num_of_plants", "rep")])
#col test
# Fit a simple model with all predictors for collinearity check
col_test <- lm(dry_root_wt ~ treatment + genotype
             + sprout_or_transplant_to_harvest, data = ds)

performance::check_collinearity(col_test)


cor(ds[, c("total_biomass", "dry_root_wt", "sprout_or_transplant_to_harvest", "num_of_plants", "rep")], use = "pairwise.complete.obs")    
pairs(ds[, c("total_biomass", "dry_root_wt", "sprout_or_transplant_to_harvest", "num_of_plants", "rep")], upper.panel = NULL)   
library(car)
vif(m2)

lm1 <- lm(dry_root_wt ~ treatment + genotype + sprout_or_transplant_to_harvest, data = ds)
summary(lm1)


# Load needed libraries
library(nlme)
library(performance)
library(sjPlot)
library(ggplot2)
library(ggeffects)

# 1️⃣ Fit GLS model (if not already done)
gls1 <- gls(dry_root_wt ~ treatment + genotype + sprout_or_transplant_to_harvest,
            data = ds,
            method = "REML",
            na.action = na.exclude)

# 2️⃣ Fit LME model with replicate as random effect
lme1 <- lme(dry_root_wt ~ treatment + genotype + sprout_or_transplant_to_harvest,
            random = ~ 1 | rep,
            data = ds,
            method = "REML",
            na.action = na.exclude)

# 3️⃣ Compare GLS vs LME with ML
gls1_ml <- update(gls1, method = "ML")
lme1_ml <- update(lme1, method = "ML")
anova(gls1_ml, lme1_ml)

# 4️⃣ Diagnostic check
check_model(lme1)

# 5️⃣ Likelihood-ratio tests for fixed effects
drop1(lme1_ml, test = "Chisq")

# 6️⃣ Final model (simplified to genotype only, using REML)
lme_final <- update(lme1, fixed = dry_root_wt ~ treatment + genotype + sprout_or_transplant_to_harvest,
                    random = ~ 1 | rep,
                    method = "REML")
summary(lme_final)
check_model(lme_final)

# 7️⃣ Present results
tab_model(lme_final, show.se = TRUE)


# Get predicted values from the model
lmer_final <- lmer(dry_root_wt ~ genotype + sprout_or_transplant_to_harvest + (1 | rep),
                   data = ds,
                   REML = TRUE)


lmer_final <- na.action(lmer_final, na.exclude)

pred <- ggpredict(lmer_final, terms = "sprout_or_transplant_to_harvest")
plot(pred, show_data = TRUE) +
  labs(x = "Sprout or Transplant to Harvest (days)",
       y = "Dry Root Weight (g)",
       title = "Predicted relationship between plant age and dry root weight") +
  theme_minimal()

library(MASS)
step_model <- stepAIC(lme1_ml, direction = "backward")
summary(step_model)

library(lme4)
library(MuMIn)

# Refit with ML
lmer_final_ml <- update(lmer_final, REML = FALSE)

# Run dredge for model selection
all_models <- dredge(lmer_final_ml)
# View the top models
all_models

# Top 5 models by AICc
head(all_models, 5)

# Best model------------------------------------------
best_model <- get.models(all_models, 1)[[1]]
summary(best_model)

library(nlme)
library(MuMIn)

lme_ml <- update(lme1, method = "ML")

all_models <- dredge(lme_ml, na.action = na.exclude)
all_models

# Average only models within delta < 4
avg_model <- model.avg(all_models, subset = delta < 4)
summary(avg_model)

performance(lme_final)

set.seed(123)
performance_cv(lme_final, method = "k_fold", k = 5, stack = FALSE)

library(lme4)

set.seed(123)





# Number of folds
k <- 5

# Get unique groups and randomly assign them to folds
groups <- unique(ds$rep)
fold_assignments <- sample(rep(1:k, length.out = length(groups)))
names(fold_assignments) <- groups

# Data frame to store results
cv_results <- data.frame(Fold = integer(), R2 = numeric(), RMSE = numeric())

# Loop through folds
for (fold in 1:k) {
  
  # Identify test and train groups
  test_groups <- names(fold_assignments[fold_assignments == fold])
  train_data <- ds[!ds$rep %in% test_groups, ]
  test_data  <- ds[ds$rep %in% test_groups, ]
  
  # Fit the mixed model on the training data
  model <- lmer(dry_root_wt ~ genotype + sprout_or_transplant_to_harvest + (1|rep),
                data = train_data)
  
  # Predict on the test data
  preds <- predict(model, newdata = test_data, allow.new.levels = TRUE)
  
  # Calculate R² and RMSE for this fold
  ss_res <- sum((test_data$dry_root_wt - preds)^2)
  ss_tot <- sum((test_data$dry_root_wt - mean(test_data$dry_root_wt))^2)
  r2 <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((test_data$dry_root_wt - preds)^2))
  
  # Store results
  cv_results <- rbind(cv_results, data.frame(Fold = fold, R2 = r2, RMSE = rmse))
}

# Average performance across folds
cv_summary <- data.frame(
  Mean_R2 = mean(cv_results$R2),
  SD_R2   = sd(cv_results$R2),
  Mean_RMSE = mean(cv_results$RMSE),
  SD_RMSE   = sd(cv_results$RMSE)
)

cv_results