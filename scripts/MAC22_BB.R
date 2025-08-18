
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
library(scales)

#Set up working directory
getwd()

#This is my path - change yours to whatever yours is
setwd("C:/Users/beabo/OneDrive/Documents/NAU/Sorghum/MAC-2022-Drought")

ds <- read.csv("data/MAC22_cleaned.csv") #Using the cleaned up dataset. Missing some things like dmgr. 

theme_set(theme_bw())

#------macombo plots----
vars <- c("shoot_wt", "florets", "d13c", "c", "d15n", "n", "p")

mycorrhizal_vars <- c("amf_in_dry_soil", "rlc_p", "dse_in_dry_soil", "amf_tot",
                      "dse_tot", "spores", "dse_hyphae", "fine_amf",
                      "coarse_amf", "ves_and_arb", "dse_and_am", "vesicle_or_spore",
                      "dse", "lse", "am_hyphae", "no_fungus",
                      "arb", "olpidium", "tot")

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
  "p" = "Phosphorus Content"#,
#  "cn_ratio" = "C:N Ratio",
 # "np_ratio" = "N:P Ratio"
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

ggsave("plots/drought_response_heatmap_simple.png", plot = p, width = 14, height = 8, dpi = 300, bg = "white")



#Multivariate Analysis here

# ---- Trait-syndrome / multivariate analysis from drought effect matrix ----

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)  # for nicer labels if available

# Build effect matrix: rows = genotype, cols = variables
effect_wide <- heatmap_df %>%
  select(genotype, variable, effect) %>%
  pivot_wider(names_from = variable, values_from = effect)

# Remove genotype column for numeric matrix, keep genotypes vector
genotypes <- effect_wide$genotype

# Ensure no Inf, then inspect NA counts
#effect_mat[is.infinite(effect_mat)] <- NA
cat("NA per column (traits):\n"); print(colSums(is.na(effect_wide)))
cat("NA per row (genotypes):\n"); print(rowSums(is.na(effect_wide)))

# Drop columns (traits) with >50% missing and rows (genotypes) with >50% missing
col_keep <- names(which(colSums(is.na(effect_wide)) <= 0.5 * nrow(effect_wide)))
row_keep <- which(rowSums(is.na(effect_wide)) <= 0.5 * ncol(effect_wide))
effect_clean <- effect_wide[row_keep, col_keep, drop = FALSE]
genotypes_clean <- genotypes[row_keep]

# ...existing code...
# Mean-impute remaining NAs per trait (ensure only numeric trait columns are used)
# remove any accidental genotype or non-numeric columns, coerce to numeric
if ("genotype" %in% names(effect_clean)) {
  effect_clean <- effect_clean %>% select(-genotype)
}

# coerce all remaining cols to numeric (warn if coercion introduces NAs)
effect_clean <- effect_clean %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

bad_cols <- names(effect_clean)[!sapply(effect_clean, is.numeric)]
if (length(bad_cols) > 0) {
  stop("Non-numeric columns present before PCA: ", paste(bad_cols, collapse = ", "))
}

cat("NA per column after coercion:\n"); print(colSums(is.na(effect_clean)))

# Mean-impute remaining NAs per trait (simple and transparent)
effect_imp <- as.data.frame(lapply(effect_clean, function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))     # fallback if all NA (should be rare)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))
rownames(effect_imp) <- genotypes_clean

# Run PCA (no na.action argument)
pca <- prcomp(effect_imp, center = TRUE, scale. = TRUE)


# Scree & loadings (as before)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
p_scree <- ggplot(data.frame(PC = seq_along(var_expl), Prop = var_expl), aes(PC, Prop)) +
  geom_line() + geom_point() +
  labs(x = "PC", y = "Proportion variance", title = "PCA scree (drought effect matrix)")
print(p_scree)

ggsave("plots/pca_scree_plot.png", plot = p_scree, width = 8, height = 6, dpi = 300)

# ...existing code...
loadings <- as.data.frame(pca$rotation)
loadings$variable <- rownames(loadings)

# prepare genotype scores to scale arrows so trait vectors overlay sensibly
scores_now <- as.data.frame(pca$x)

# compute multiplier to scale loadings to score space (safe to NA)
mult <- 0.7 * max(abs(c(scores_now$PC1, scores_now$PC2)), na.rm = TRUE) /
         max(abs(c(loadings$PC1, loadings$PC2)), na.rm = TRUE)

library(grid)  # for unit()

p_load <- ggplot() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  # optional: show genotype scores in the background for context
  geom_point(data = scores_now, aes(x = PC1, y = PC2), color = "gray60", size = 2, alpha = 0.6) +
  # arrows from origin to scaled loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.18, "cm")), color = "firebrick", alpha = 0.9) +
  # labels at arrow tips
  ggrepel::geom_text_repel(data = loadings,
                           aes(x = PC1 * mult, y = PC2 * mult, label = variable),
                           color = "firebrick", size = 3) +
  coord_equal() +
  labs(title = "Trait loadings on PC1 vs PC2 (drought effects)", x = "PC1", y = "PC2") +
  theme_minimal()

print(p_load)


ggsave("plots/pca_trait_loadings.png", plot = p_load, width = 8, height = 6, dpi = 300)

# ...existing code...
# compute simple sensitivity (mean effect across traits) and add to scores_now
sensitivity <- rowMeans(effect_imp, na.rm = TRUE)   # effect_imp rows = genotypes
scores_now$sensitivity <- sensitivity[rownames(scores_now)]

p_load_colored <- p_load + 
  geom_point(data = scores_now, aes(x = PC1, y = PC2, color = sensitivity), size = 3, inherit.aes = FALSE) +
  scale_color_gradient2(midpoint = 0, low = "#2166AC", mid = "white", high = "#B2182B", name = "Mean effect")
print(p_load_colored)
ggsave("plots/pca_trait_loadings_colored_by_sensitivity.png", plot = p_load_colored, width = 8, height = 6, dpi = 300)


# build genotype × treatment matrix (mean per cell), then PCA
gt_mat <- filtered_long_data %>%
  group_by(genotype, treatment, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variable, values_from = mean_val)

# keep ids and numeric matrix
gt_ids <- gt_mat %>% select(genotype, treatment)
gt_num <- gt_mat %>% select(-genotype, -treatment) %>% mutate(across(everything(), ~ as.numeric(.)))
# simple mean impute if needed
gt_num <- as.data.frame(lapply(gt_num, function(x){ x[is.na(x)] <- mean(x, na.rm=TRUE); x }))
rownames(gt_num) <- paste(gt_ids$genotype, gt_ids$treatment, sep = "_")

pca_gt <- prcomp(gt_num, center = TRUE, scale. = TRUE)
scores_gt <- as.data.frame(pca_gt$x)
scores_gt$id <- rownames(scores_gt)
scores_gt <- tidyr::separate(scores_gt, id, into = c("genotype", "treatment"), sep = "_", remove = FALSE)

# plot PC1 vs PC2 colored by treatment
p_gt <- ggplot(scores_gt, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  labs(title = "PCA of genotype × treatment profiles", color = "Treatment") +
  theme_minimal()
print(p_gt)
ggsave("plots/pca_genotype_treatment_by_treatment.png", plot = p_gt, width = 8, height = 6, dpi = 300)

# assume pca_gt exists and scores_gt was created
loadings_gt <- as.data.frame(pca_gt$rotation)
loadings_gt$trait <- rownames(loadings_gt)

scores_gt <- as.data.frame(pca_gt$x)
scores_gt$id <- rownames(scores_gt)
scores_gt <- tidyr::separate(scores_gt, id, into = c("genotype","treatment"), sep = "_", remove = FALSE)

# scale loadings into score space
mult <- 0.7 * max(abs(c(scores_gt$PC1, scores_gt$PC2)), na.rm = TRUE) /
        max(abs(c(loadings_gt$PC1, loadings_gt$PC2)), na.rm = TRUE)

library(ggrepel); library(grid)
p_gt_biplot <- ggplot() +
  geom_point(data = scores_gt, aes(PC1, PC2, color = treatment), size = 3, alpha = 0.8) +
  geom_segment(data = loadings_gt,
               aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult),
               arrow = arrow(length = unit(0.18, "cm")), color = "darkred") +
  geom_text_repel(data = loadings_gt, 
                  aes(x = PC1 * mult, y = PC2 * mult, label = trait),
                  color = "darkred", size = 3) +
  coord_equal() + theme_minimal() +
  labs(title = "PCA (genotype:treatment) with trait loadings")
print(p_gt_biplot)

ggsave("plots/pca_genotype_treatment_biplot.png", plot = p_gt_biplot, width = 8, height = 6, dpi = 300)

# Build genotype scores and a performance dataframe using trait-specific effects (shoot_wt, florets)
scores <- as.data.frame(pca$x)
scores$genotype <- rownames(scores)

perf <- data.frame(genotype = rownames(effect_imp))
if ("shoot_wt" %in% colnames(effect_imp)) perf$shoot_wt_eff <- effect_imp[, "shoot_wt"] else perf$shoot_wt_eff <- NA
if ("florets" %in% colnames(effect_imp)) perf$florets_eff <- effect_imp[, "florets"] else perf$florets_eff <- NA
# Add P effect if present in the effect matrix
if ("p" %in% colnames(effect_imp)) perf$p_effect <- effect_imp[, "p"] else perf$p_effect <- NA

scores_perf <- dplyr::left_join(scores, perf, by = "genotype")

# Fit simple linear models: performance ~ PC1 + PC2
if (!all(is.na(scores_perf$shoot_wt_eff))) {
  lm_shoot <- lm(shoot_wt_eff ~ PC1 + PC2, data = scores_perf)
  cat("Shoot weight model:\n"); print(summary(lm_shoot))
  p_shoot <- ggplot(scores_perf, aes(x = PC1, y = shoot_wt_eff, label = genotype)) +
    geom_point() + geom_smooth(method = "lm", se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs shoot_wt effect")
  print(p_shoot)
}


ggsave("plots/pca_shoot_weight_effect.png", p_shoot, width = 8, height = 6, dpi = 300)

if (!all(is.na(scores_perf$florets_eff))) {
  lm_florets <- lm(florets_eff ~ PC1 + PC2, data = scores_perf)
  cat("Florets model:\n"); print(summary(lm_florets))
  p_florets <- ggplot(scores_perf, aes(x = PC1, y = florets_eff, label = genotype)) +
    geom_point() + geom_smooth(method = "lm", se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs florets effect")
  print(p_florets)
}

ggsave("plots/pca_florets_effect.png", p_florets, width = 8, height = 6, dpi = 300)

# p (phosphorus) effect: model and plot if present
if (!all(is.na(scores_perf$p_effect))) {
  lm_p <- lm(p_effect ~ PC1 + PC2, data = scores_perf)
  cat("P effect model:\n"); print(summary(lm_p))
  p_p <- ggplot(scores_perf, aes(x = PC1, y = p_effect, label = genotype)) +
    geom_point() + geom_smooth(method = "lm", se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs P effect")
  print(p_p)
  ggsave("plots/pca_p_effect.png", p_p, width = 8, height = 6, dpi = 300)
}

# ---- mycorrhizal variables correlation with PC1 ----
mycorrhizal_vars <- c("amf_in_dry_soil", "rlc_p", "dse_in_dry_soil", "amf_tot",
                      "dse_tot", "spores", "dse_hyphae", "fine_amf",
                      "coarse_amf", "ves_and_arb", "dse_and_am", "vesicle_or_spore",
                      "dse", "lse", "am_hyphae", "no_fungus",
                      "arb", "olpidium", "tot")

# Compute genotype-level drought effect for each mycorrhizal var (Droughted - Watered)
myco_long <- ds %>% select(genotype, treatment, intersect(mycorrhizal_vars, names(ds))) %>%
  pivot_longer(cols = -c(genotype, treatment), names_to = "variable", values_to = "value") %>%
  group_by(genotype, treatment, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = treatment, values_from = mean_val) %>%
  mutate(effect = Droughted - Watered) %>%
  filter(!is.na(effect))

if (nrow(myco_long) > 0) {
  myco_wide <- myco_long %>% select(genotype, variable, effect) %>% pivot_wider(names_from = variable, values_from = effect)
  # align with PCA genotypes
  common_genos <- intersect(rownames(effect_imp), myco_wide$genotype)
  if (length(common_genos) > 0) {
  myco_sub <- myco_wide %>% filter(genotype %in% common_genos)
  # convert to data.frame (avoid tibble rowname deprecation) and set rownames
  myco_sub <- as.data.frame(myco_sub)
  rownames(myco_sub) <- myco_sub$genotype
  myco_sub$genotype <- NULL
  # get PC1 scores matched to genotype order
  pc1_scores <- scores$PC1[match(rownames(myco_sub), scores$genotype)]
    cor_list <- lapply(colnames(myco_sub), function(var) {
      eff <- as.numeric(myco_sub[[var]])
      ok <- !is.na(eff) & !is.na(pc1_scores)
      if (sum(ok) >= 3) {
        ct <- cor.test(eff[ok], pc1_scores[ok])
        data.frame(variable = var, cor = unname(ct$estimate), p.value = ct$p.value, n = sum(ok))
      } else {
        data.frame(variable = var, cor = NA_real_, p.value = NA_real_, n = sum(ok))
      }
    })
    cor_df <- do.call(rbind, cor_list)
    cor_df$padj <- p.adjust(cor_df$p.value, method = 'BH')
    cor_df_ord <- cor_df %>% arrange(desc(abs(cor)))
    write.csv(cor_df_ord, file.path('results', 'pca_mycorrhizal_correlations_with_fdr.csv'), row.names = FALSE)

    # Focus on user-selected mycorrhizal variables (fallback to available significant vars if none present)
    preferred_vars <- c("dse_tot", "amf_tot")
    avail_pref <- intersect(preferred_vars, colnames(myco_sub))

    if (length(avail_pref) > 0) {
      sel_vars <- avail_pref
      cat('Plotting user-selected mycorrhizal vars:', paste(sel_vars, collapse = ', '), '\n')
    } else {
      # fallback: use top-correlated vars (by absolute correlation)
      sel_vars <- head(cor_df_ord$variable, 6)
      cat('No preferred mycorrhizal vars present; falling back to top correlated vars:', paste(sel_vars, collapse = ', '), '\n')
    }

    if (length(sel_vars) == 0) {
      cat('No mycorrhizal variables available for plotting. Skipping combined plot.\n')
    } else {
      combined_df <- do.call(rbind, lapply(sel_vars, function(var) {
        data.frame(genotype = rownames(myco_sub), effect = as.numeric(myco_sub[[var]]), PC1 = pc1_scores, variable = var, stringsAsFactors = FALSE)
      }))
      combined_df <- combined_df[!is.na(combined_df$effect) & !is.na(combined_df$PC1), ]

      p_comb <- ggplot(combined_df, aes(x = PC1, y = effect, color = variable)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_smooth(method = 'lm', se = FALSE, linetype = 'solid', size = 1) +
        labs(title = 'PC1 vs selected mycorrhizal effects', x = 'PC1', y = 'Effect (Drought - Watered)') +
        theme_minimal() +
        theme(legend.title = element_blank())

      print(p_comb)
      ggsave(file.path('plots', 'pca_myco_selected_vs_PC1.png'), plot = p_comb, width = 10, height = 6, dpi = 300)
    }
  }
}


# Save outputs
write.csv(loadings, "results/pca_trait_loadings.csv", row.names = FALSE)
write.csv(scores_perf, "results/pca_genotype_scores_and_performance.csv", row.names = FALSE)


#Correlate mycorrhizal traits to PC1.



#STOP HERE

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
