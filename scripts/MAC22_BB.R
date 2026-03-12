
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

# Add group ellipses based on trait_group_k / cluster_for_geno (if available)
try({
  # ensure scores_now has genotype identifiers matching effect_imp / scores
  if (!"genotype" %in% colnames(scores_now)) scores_now$genotype <- rownames(scores_now)
  # attach group labels where available
  if (exists("cluster_for_geno") && length(cluster_for_geno) == nrow(effect_imp)) {
    gmap <- data.frame(genotype = rownames(effect_imp), group = as.factor(cluster_for_geno), stringsAsFactors = FALSE)
    scores_now <- merge(scores_now, gmap, by = 'genotype', all.x = TRUE)
  } else if (exists("scores") && "genotype" %in% colnames(scores)) {
    scores_now <- merge(scores_now, scores[, c('genotype')], by.x = 'row.names', by.y = 'genotype', all.x = TRUE)
    rownames(scores_now) <- scores_now$Row.names; scores_now$Row.names <- NULL
    if (exists("cluster_for_geno")) scores_now$group <- as.factor(cluster_for_geno[match(rownames(scores_now), names(cluster_for_geno))])
  }

  if ("group" %in% colnames(scores_now)) {
    # draw ellipse fill and border, slightly transparent
    p_ell <- p_load +
      geom_point(data = scores_now, aes(x = PC1, y = PC2, color = group), size = 2, inherit.aes = FALSE) +
      stat_ellipse(data = scores_now, aes(x = PC1, y = PC2, fill = group, color = group), geom = 'polygon', alpha = 0.15, level = 0.68, show.legend = TRUE) +
      scale_fill_brewer(palette = 'Set2') +
      scale_color_brewer(palette = 'Set2')
    ggsave("plots/pca_trait_loadings_with_group_ellipses.png", plot = p_ell, width = 8, height = 6, dpi = 300)
  }
}, silent = TRUE)

# Also produce a version with exactly two group ellipses (visual only)
try({
  s2 <- scores_now
  if (!"genotype" %in% colnames(s2)) s2$genotype <- rownames(s2)
  # if cluster_for_geno exists and has length matching effect_imp, prefer mapping; otherwise recluster PC1 into 2 groups
  if (exists("cluster_for_geno") && length(cluster_for_geno) == nrow(effect_imp)) {
    # create binary grouping by kmeans on PC1 to avoid altering saved cluster_for_geno
    km2 <- kmeans(s2$PC1, centers = 2)
    s2$group2 <- as.factor(km2$cluster)
  } else {
    km2 <- kmeans(s2$PC1, centers = 2)
    s2$group2 <- as.factor(km2$cluster)
  }

  p_2g <- p_load +
    geom_point(data = s2, aes(x = PC1, y = PC2, color = group2), size = 2, inherit.aes = FALSE) +
    stat_ellipse(data = s2, aes(x = PC1, y = PC2, fill = group2, color = group2), geom = 'polygon', alpha = 0.15, level = 0.68, show.legend = TRUE) +
    scale_fill_brewer(palette = 'Set1') +
    scale_color_brewer(palette = 'Set1') +
    labs(title = 'Trait loadings on PC1 vs PC2 (2-group ellipses)')
  ggsave('plots/pca_trait_loadings_with_2group_ellipses.png', plot = p_2g, width = 8, height = 6, dpi = 300)
}, silent = TRUE)


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
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs shoot_wt effect")
  print(p_shoot)
}


ggsave("plots/pca_shoot_weight_effect.png", p_shoot, width = 8, height = 6, dpi = 300)

if (!all(is.na(scores_perf$florets_eff))) {
  lm_florets <- lm(florets_eff ~ PC1 + PC2, data = scores_perf)
  cat("Florets model:\n"); print(summary(lm_florets))
  p_florets <- ggplot(scores_perf, aes(x = PC1, y = florets_eff, label = genotype)) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs florets effect")
  print(p_florets)
}

ggsave("plots/pca_florets_effect.png", p_florets, width = 8, height = 6, dpi = 300)

# p (phosphorus) effect: model and plot if present
if (!all(is.na(scores_perf$p_effect))) {
  lm_p <- lm(p_effect ~ PC1 + PC2, data = scores_perf)
  cat("P effect model:\n"); print(summary(lm_p))
  p_p <- ggplot(scores_perf, aes(x = PC1, y = p_effect, label = genotype)) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + ggrepel::geom_text_repel(size = 3) +
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

# ---- PC1 backing metrics ----
# Variance explained (PC1)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_var_pct <- 100 * var_expl[1]

# Top trait loadings on PC1
top_loadings <- data.frame(variable = rownames(pca$rotation), loading = pca$rotation[,"PC1"], stringsAsFactors = FALSE)
top_loadings <- top_loadings[order(-abs(top_loadings$loading)), ]
write.csv(top_loadings, "results/pca_PC1_top_loadings.csv", row.names = FALSE)

# Define PC1 groups (tertiles) and per-group summaries
scores_perf$PC1_group <- dplyr::ntile(scores_perf$PC1, 2)
group_summary <- scores_perf %>% dplyr::group_by(PC1_group) %>% dplyr::summarise(n = dplyr::n(), mean = mean(PC1, na.rm=TRUE), sd = sd(PC1, na.rm=TRUE), se = sd/sqrt(n), lo = mean - 1.96*se, hi = mean + 1.96*se)
write.csv(group_summary, "results/pca_PC1_group_summary.csv", row.names = FALSE)

# ANOVA / linear model: PC1 ~ group
lm_group <- lm(PC1 ~ factor(PC1_group), data = scores_perf)
anova_group <- anova(lm_group)
group_r2 <- summary(lm_group)$r.squared
write.csv(data.frame(stat = c('R2'), value = c(group_r2)), "results/pca_PC1_group_r2.csv", row.names = FALSE)

# Cohen's d for top vs bottom tertiles (if present)
cohens_d <- NA_real_
if (all(c(1,2) %in% scores_perf$PC1_group)) {
  g1 <- scores_perf$PC1[scores_perf$PC1_group == 1]
  g2 <- scores_perf$PC1[scores_perf$PC1_group == 2]
  pooled_sd <- sqrt(((length(g1)-1)*var(g1, na.rm=TRUE) + (length(g2)-1)*var(g2, na.rm=TRUE)) / (length(g1)+length(g2)-2))
  cohens_d <- (mean(g2, na.rm=TRUE) - mean(g1, na.rm=TRUE)) / pooled_sd
}
write.csv(data.frame(cohens_d = cohens_d), "results/pca_PC1_cohens_d_top_bottom.csv", row.names = FALSE)

# Permutation test for top vs bottom tertiles (difference in means)
perm_p <- NA_real_
if (all(c(1,3) %in% scores_perf$PC1_group)) {
  obs_diff <- mean(scores_perf$PC1[scores_perf$PC1_group==3], na.rm=TRUE) - mean(scores_perf$PC1[scores_perf$PC1_group==1], na.rm=TRUE)
  perm_diffs <- replicate(2000, {
    perm_group <- sample(scores_perf$PC1_group)
    mean(scores_perf$PC1[perm_group==3], na.rm=TRUE) - mean(scores_perf$PC1[perm_group==1], na.rm=TRUE)
  })
  perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))
}
write.csv(data.frame(perm_p = perm_p), "results/pca_PC1_perm_p_top_bottom.csv", row.names = FALSE)

# Bootstrap CI for group means (example for group 1 and 3)
boot_ci <- function(x, B = 2000) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(c(NA, NA))
  sims <- replicate(B, mean(sample(x, replace = TRUE)))
  quantile(sims, c(0.025, 0.975))
}
ci1 <- boot_ci(scores_perf$PC1[scores_perf$PC1_group==1])
ci2 <- boot_ci(scores_perf$PC1[scores_perf$PC1_group==2])
write.csv(data.frame(group = c(1,2), lo = c(ci1[1], ci2[1]), hi = c(ci1[2], ci2[2])), "results/pca_PC1_bootstrap_CIs_top_bottom.csv", row.names = FALSE)

# PERMANOVA across the multivariate effect matrix if vegan available
if (requireNamespace('vegan', quietly = TRUE)) {
  groups_for_adonis <- scores_perf$PC1_group[match(rownames(effect_imp), scores_perf$genotype)]
  if (length(groups_for_adonis) == nrow(effect_imp)) {
    ad <- vegan::adonis2(as.matrix(effect_imp) ~ groups_for_adonis, permutations = 999)
    capture_out <- capture.output(print(ad))
    writeLines(capture_out, con = "results/pca_PC1_permanova.txt")
  }
}

# Correlations: PC1 with selected variables (shoot, florets, p, and some mycorrhizal vars)
vars_to_corr <- c('shoot_wt_eff', 'florets_eff', 'p_effect', 'dse_hyphae', 'dse_in_dry_soil', 'dse_tot', 'amf_tot', 'lse', 'no_fungus')
corr_res <- lapply(vars_to_corr, function(v) {
  y <- scores_perf[[v]]
  ok <- !is.na(scores$PC1) & !is.na(y[match(scores$genotype, scores_perf$genotype)])
  x <- scores$PC1[ok]
  yy <- y[match(scores$genotype, scores_perf$genotype)][ok]
  if (length(x) >= 3) {
    ct <- cor.test(x, yy)
    data.frame(variable = v, cor = as.numeric(ct$estimate), p.value = ct$p.value, n = length(x))
  } else data.frame(variable = v, cor = NA_real_, p.value = NA_real_, n = length(x))
})
corr_df <- do.call(rbind, corr_res)
corr_df$padj <- p.adjust(corr_df$p.value, method = 'BH')
write.csv(corr_df, 'results/pca_PC1_correlations.csv', row.names = FALSE)

# Append a short numeric summary to the plot_summaries log
sink('results/plot_summaries.txt', append = TRUE)
cat('\n---- PC1 backing metrics ----\n')
cat(sprintf('PC1 variance (%%): %.2f\n', pc1_var_pct))
cat('Top loadings (top 6):\n'); print(head(top_loadings, 6))
cat('\nGroup summary (PC1 tertiles):\n'); print(group_summary)
cat('\nANOVA (PC1 ~ tertile):\n'); print(anova_group)
cat('\nR-squared (PC1 ~ tertile):\n'); print(group_r2)
cat('\nCohen\'s d (top vs bottom tertiles):\n'); print(cohens_d)
cat('\nPermutation p (top vs bottom tertiles):\n'); print(perm_p)
cat('\nCorrelation summary (PC1 vs variables):\n'); print(corr_df)
sink()

# ---- Trait grouping by PC1 groups (automatic k selection) ----
# Choose k automatically (2..max_k) by maximizing average silhouette width on PC1
set.seed(1998)
# restrict maximum tested k to 3 by default to avoid over-splitting on small datasets
max_k <- min(3, max(2, nrow(scores) - 1))
chosen_k <- 2
if (nrow(scores) >= 2) {
  possible_k <- seq(2, max_k)
  sil_scores <- rep(NA_real_, length(possible_k))
  if (requireNamespace('cluster', quietly = TRUE)) {
    dmat <- dist(as.matrix(scores$PC1))
    min_cluster_size <- 2
    sil_log <- list()
    for (i in seq_along(possible_k)) {
      k <- possible_k[i]
      km <- tryCatch(kmeans(scores$PC1, centers = k, nstart = 25), error = function(e) NULL)
      if (is.null(km)) next
      labs <- km$cluster
      cluster_sizes <- as.integer(table(labs))
      # require minimum cluster size to consider this k
      if (any(cluster_sizes < min_cluster_size)) {
        sil_scores[i] <- NA_real_
        sil_log[[as.character(k)]] <- list(sil = NA_real_, sizes = cluster_sizes)
        next
      }
      if (length(unique(labs)) < 2) next
      sil <- tryCatch(cluster::silhouette(labs, dmat), error = function(e) NULL)
      if (is.null(sil)) next
      sil_scores[i] <- mean(sil[, 3], na.rm = TRUE)
      sil_log[[as.character(k)]] <- list(sil = sil_scores[i], sizes = cluster_sizes)
    }
    # write silhouette log
    try({
      sink(file.path('results', 'pca_silhouette_log.txt'), append = TRUE)
      cat('silhouette log for possible ks:\n')
      print(sil_log)
      sink()
    }, silent = TRUE)
    if (all(is.na(sil_scores))) {
      chosen_k <- 2
    } else {
      chosen_k <- possible_k[which.max(sil_scores)]
    }
  } else {
    # fallback when cluster not available
    chosen_k <- 2
  }
} else {
  chosen_k <- 1
}
# Save chosen k
write.csv(data.frame(chosen_k = chosen_k), file = file.path('results', 'pca_chosen_k.csv'), row.names = FALSE)
# use chosen_k for downstream grouping
trait_group_k <- chosen_k
if (nrow(scores) >= trait_group_k) {
  if (trait_group_k == 1) {
    # single cluster: assign all genotypes to group 1
    raw_clusters <- rep(1, nrow(scores))
    cm <- mean(scores$PC1, na.rm = TRUE)
    clusters_ord <- raw_clusters
  } else {
    raw_clusters <- kmeans(scores$PC1, centers = trait_group_k)$cluster
    # reorder cluster labels by ascending cluster mean for interpretability
  cm <- tapply(scores$PC1, raw_clusters, mean)
  ord <- order(cm)
  # remap original cluster labels to ordered ranks (smallest mean -> 1)
  remap <- setNames(seq_along(ord), names(cm)[ord])
  clusters_ord <- remap[as.character(raw_clusters)]
  }

  # align clusters to effect_imp genotypes
  cluster_for_geno <- clusters_ord[match(rownames(effect_imp), scores$genotype)]

  # robustly compact labels so observed groups become consecutive 1..m (avoid gaps like 1 and 3)
  # ensure cluster_for_geno is an integer vector (may contain NA)
  cluster_for_geno <- as.integer(cluster_for_geno)
  observed <- sort(unique(cluster_for_geno[!is.na(cluster_for_geno)]))
  if (length(observed) == 0) {
    # no clusters observed: mark all as NA and treat as single-group downstream
    cluster_for_geno[] <- NA_integer_
    trait_group_k <- 1
  } else {
    # build mapping from original observed labels -> consecutive 1..m
    map_obs <- setNames(seq_along(observed), as.character(observed))
    # map values safely, preserving NAs
    cluster_chr <- as.character(cluster_for_geno)
    mapped <- map_obs[cluster_chr]
    # any entries not found in map_obs will become NA; coerce to integer
    cluster_for_geno <- as.integer(mapped)
    # update trait_group_k to the actual number of observed groups after mapping
    trait_group_k <- length(observed)
  }

  # Diagnostic logging: help debug why group means may be all NA
  dbg_file <- file.path('results', paste0('pca_trait_group_debug_k', trait_group_k, '.txt'))
  dir.create(dirname(dbg_file), showWarnings = FALSE, recursive = TRUE)
  sink(dbg_file)
  cat('Diagnostic for trait grouping by PC1\n')
  cat('trait_group_k =', trait_group_k, '\n')
  cat('nrow(effect_imp) =', nrow(effect_imp), '\n')
  cat('length(scores$genotype) =', length(scores$genotype), '\n')
  cat('first 10 effect_imp rownames:\n'); print(head(rownames(effect_imp), 10))
  cat('first 10 scores genotypes:\n'); print(head(scores$genotype, 10))
  matched <- rownames(effect_imp) %in% scores$genotype
  cat('matched count:', sum(matched), 'of', length(matched), '\n')
  if (any(!matched)) {
    cat('unmatched (first 20):\n'); print(head(rownames(effect_imp)[!matched], 20))
  }
  cat('raw_clusters length:', length(raw_clusters), '\n')
  cat('clusters_ord length:', length(clusters_ord), '\n')
  cat('cluster_for_geno length:', length(cluster_for_geno), '\n')
  cat('cluster_for_geno summary:\n'); print(summary(cluster_for_geno))
  sink()

  trait_stats_list <- lapply(colnames(effect_imp), function(tr) {
    vals <- effect_imp[[tr]]
    # ensure we only consider genotypes with a cluster assignment
    ok_idx <- !is.na(cluster_for_geno)
    groups_present <- unique(cluster_for_geno[ok_idx & !is.na(vals)])
    # initialize full-length group means vector
    group_means_full <- rep(NA_real_, trait_group_k)
    names(group_means_full) <- seq_len(trait_group_k)
    if (length(groups_present) > 0) {
      mm <- tapply(vals[ok_idx], cluster_for_geno[ok_idx], mean, na.rm = TRUE)
      # mm may be named by group numbers; fill into full vector
      if (length(mm) > 0) {
        group_means_full[names(mm)] <- as.numeric(mm)
      }
    }
    group_ns <- sapply(seq_len(trait_group_k), function(g) sum(!is.na(vals) & cluster_for_geno == g))
    pval <- NA_real_
    if (trait_group_k == 2) {
      # two-group t-test when both groups have >=2 non-missing
      if (all(group_ns >= 2, na.rm = TRUE)) {
        g1 <- vals[cluster_for_geno == 1]
        g2 <- vals[cluster_for_geno == 2]
        if (sum(!is.na(g1)) >= 2 && sum(!is.na(g2)) >= 2) {
          pval <- tryCatch(t.test(g1, g2)$p.value, error = function(e) NA_real_)
        }
      }
    } else if (trait_group_k > 2) {
      # For >2 groups, try ANOVA (requires at least two groups with >=2 samples); fallback to Kruskal-Wallis
      df_tv <- data.frame(val = vals, grp = cluster_for_geno)
      df_tv <- df_tv[!is.na(df_tv$val) & !is.na(df_tv$grp), , drop = FALSE]
      if (nrow(df_tv) >= 3) {
        gs <- table(df_tv$grp)
        if (sum(gs >= 2) >= 2) {
          # try ANOVA
          a_mod <- tryCatch(aov(val ~ factor(grp), data = df_tv), error = function(e) NULL)
          if (!is.null(a_mod)) {
            an <- tryCatch(anova(a_mod), error = function(e) NULL)
            if (!is.null(an) && nrow(an) >= 1) {
              pval <- as.numeric(an$"Pr(>F)"[1])
            }
          }
          # fallback to Kruskal-Wallis if ANOVA failed or returned NA
          if (is.na(pval)) {
            kw <- tryCatch(kruskal.test(val ~ factor(grp), data = df_tv), error = function(e) NULL)
            if (!is.null(kw)) pval <- kw$p.value
          }
        }
      }
    }
    if (all(is.na(group_means_full))) assigned <- NA_integer_ else assigned <- as.integer(which.max(group_means_full))
    data.frame(trait = tr, assigned_group = assigned, p.value = pval, stringsAsFactors = FALSE)
  })
  trait_stats <- do.call(rbind, trait_stats_list)
  trait_stats$padj <- p.adjust(trait_stats$p.value, method = 'BH')
  write.csv(trait_stats, file.path('results', paste0('pca_trait_group_stats_PC1_k', trait_group_k, '.csv')), row.names = FALSE)

  # write per-trait group means table
  means_mat <- t(sapply(colnames(effect_imp), function(tr) {
    mm <- tapply(effect_imp[[tr]], cluster_for_geno, mean, na.rm = TRUE)
    # ensure length equals trait_group_k
    out <- rep(NA_real_, trait_group_k)
    if (length(mm) > 0) {
      # names(mm) are group numbers like '1','2' — assign by numeric index
      idx <- as.integer(names(mm))
      out[idx] <- as.numeric(mm)
    }
    names(out) <- paste0('group', seq_len(trait_group_k))
    out
  }))
  means_df <- data.frame(trait = rownames(means_mat), means_mat, stringsAsFactors = FALSE, row.names = NULL)
  write.csv(means_df, file.path('results', paste0('pca_trait_group_means_PC1_k', trait_group_k, '.csv')), row.names = FALSE)

  # simple heatmap of trait means by group
  long_means <- tidyr::pivot_longer(means_df, cols = starts_with('group'), names_to = 'group', values_to = 'mean')
  p_heat <- ggplot(long_means, aes(x = group, y = trait, fill = mean)) + geom_tile() + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', na.value = 'grey50') + labs(title = paste('Trait means by PC1 groups (k=', trait_group_k, ')', sep='')) + theme_minimal()
  ggsave(file.path('plots', paste0('pca_PC1_trait_group_means_k', trait_group_k, '.png')), plot = p_heat, width = 8, height = max(6, 0.25 * nrow(means_df)), dpi = 300)

  # append short summary
  sink('results/plot_summaries.txt', append = TRUE)
  cat('\n---- Trait grouping by PC1 (k=', trait_group_k, ') ----\n', sep='')
  cat('Trait group counts:\n')
  print(table(trait_stats$assigned_group))
  cat('\nTop traits assigned to each group (by mean):\n')
  for (g in seq_len(trait_group_k)) {
    topg <- means_df %>% dplyr::arrange(dplyr::desc(.data[[paste0('group', g)]])) %>% dplyr::slice_head(n = 8) %>% dplyr::pull(trait)
    cat('Group', g, ':', paste(topg, collapse = ', '), '\n')
  }
  sink()
}

#Correlate mycorrhizal traits to PC1.

stop("stop")

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
