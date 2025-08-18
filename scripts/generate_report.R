# generate_report.R
# Recompute key summaries and build a markdown report that embeds plots and numeric outputs

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Paths
data_file <- "data/MAC22_cleaned.csv"
plots_dir <- "plots"
results_dir <- "results"
summary_txt <- file.path(results_dir, "plot_summaries.txt")
report_md <- file.path(results_dir, "report.md")

# Ensure results dir
if(!dir.exists(results_dir)) dir.create(results_dir)

# Read cleaned data
ds <- read.csv(data_file, stringsAsFactors = FALSE)

vars <- c("shoot_wt", "florets", "d13c", "c", "d15n", "n", "p")

# scale by genotype (same as main script)
ds_scaled <- ds %>% group_by(genotype) %>% mutate(across(all_of(vars), ~ scale(.)[,1])) %>% ungroup()

# long format
filtered_long_data <- ds_scaled %>% select(genotype, treatment, all_of(vars)) %>% pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>% filter(!is.na(value))

# safe t-tests per genotype-variable
stats_df <- filtered_long_data %>% group_by(genotype, variable) %>% summarise(p_value = tryCatch({ if(length(unique(treatment))==2) t.test(value ~ treatment)$p.value else NA_real_ }, error = function(e) NA_real_), .groups = 'drop') %>% mutate(stars = case_when(!is.na(p_value) & p_value <= 0.001 ~ '***', !is.na(p_value) & p_value <= 0.01 ~ '**', !is.na(p_value) & p_value <= 0.05 ~ '*', TRUE ~ ''))

# heatmap effect (Droughted - Watered)
heatmap_df <- filtered_long_data %>% group_by(genotype, treatment, variable) %>% summarise(mean_val = mean(value, na.rm=TRUE), .groups='drop') %>% pivot_wider(names_from = treatment, values_from = mean_val) %>% mutate(effect = Droughted - Watered) %>% filter(!is.na(effect))

# Basic numeric summaries
n_significant <- sum(stats_df$p_value < 0.05, na.rm = TRUE)
top_effects <- heatmap_df %>% mutate(abs_eff = abs(effect)) %>% arrange(desc(abs_eff)) %>% slice_head(n = 10) %>% select(genotype, variable, effect)

# Rebuild effect matrix for PCA (same logic as main script)
effect_wide <- heatmap_df %>% select(genotype, variable, effect) %>% pivot_wider(names_from = variable, values_from = effect)

genotypes <- effect_wide$genotype

# handle missing
col_keep <- names(which(colSums(is.na(effect_wide)) <= 0.5 * nrow(effect_wide)))
row_keep <- which(rowSums(is.na(effect_wide)) <= 0.5 * ncol(effect_wide))
effect_clean <- effect_wide[row_keep, col_keep, drop = FALSE]
genotypes_clean <- genotypes[row_keep]

if ("genotype" %in% names(effect_clean)) effect_clean <- effect_clean %>% select(-genotype)

effect_clean <- effect_clean %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

effect_imp <- as.data.frame(lapply(effect_clean, function(x){ if(all(is.na(x))) return(rep(0, length(x))); x[is.na(x)] <- mean(x, na.rm = TRUE); x }))
rownames(effect_imp) <- genotypes_clean

# PCA
pca <- prcomp(effect_imp, center = TRUE, scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
loadings <- as.data.frame(pca$rotation)
loadings$variable <- rownames(loadings)

# Save numeric PCA outputs
write.csv(loadings, file.path(results_dir, "pca_trait_loadings_regen.csv"), row.names = FALSE)
write.csv(data.frame(PC = seq_along(var_expl), VarExplained = var_expl, Cumulative = cumsum(var_expl)), file.path(results_dir, "pca_variance_explained.csv"), row.names = FALSE)

# Build scores and performance
scores <- as.data.frame(pca$x)
scores$genotype <- rownames(scores)
perf <- data.frame(genotype = rownames(effect_imp))
perf$shoot_wt_eff <- if("shoot_wt" %in% colnames(effect_imp)) effect_imp[, "shoot_wt"] else NA
perf$florets_eff <- if("florets" %in% colnames(effect_imp)) effect_imp[, "florets"] else NA
perf$p_effect <- if("p" %in% colnames(effect_imp)) effect_imp[, "p"] else NA
scores_perf <- left_join(scores, perf, by = "genotype")
write.csv(scores_perf, file.path(results_dir, "pca_genotype_scores_and_performance_regen.csv"), row.names = FALSE)

# Linear models
lm_shoot <- NULL; lm_florets <- NULL; lm_p <- NULL; lm_shoot_sum <- NULL; lm_florets_sum <- NULL; lm_p_sum <- NULL
if (!all(is.na(scores_perf$shoot_wt_eff))) { lm_shoot <- lm(shoot_wt_eff ~ PC1 + PC2, data = scores_perf); lm_shoot_sum <- summary(lm_shoot) }
if (!all(is.na(scores_perf$florets_eff))) { lm_florets <- lm(florets_eff ~ PC1 + PC2, data = scores_perf); lm_florets_sum <- summary(lm_florets) }
if (!all(is.na(scores_perf$p_effect))) { lm_p <- lm(p_effect ~ PC1 + PC2, data = scores_perf); lm_p_sum <- summary(lm_p) }

# Optional: create a simple scatter+lm plot for p_effect
if (!all(is.na(scores_perf$p_effect))) {
  p_p <- ggplot(scores_perf, aes(x = PC1, y = p_effect, label = genotype)) +
    geom_point() + geom_smooth(method = "lm", se = TRUE) + ggrepel::geom_text_repel(size = 3) +
    labs(title = "PC1 vs P effect")
  ggsave(file.path(plots_dir, "pca_p_effect.png"), plot = p_p, width = 8, height = 6, dpi = 300)
}

# Append numeric outputs to results/plot_summaries.txt
sink(summary_txt, append = TRUE)
cat('\n\n--- Numeric summaries appended by scripts/generate_report.R ---\n')
cat('Number of significant genotype x trait tests (p < 0.05): ', n_significant, '\n')
cat('\nTop 10 absolute effects (genotype, variable, effect):\n')
print(top_effects)
cat('\nPCA variance explained (first 6 PCs):\n')
print(round(head(data.frame(PC = seq_along(var_expl), VarPct = 100*var_expl, CumPct = 100*cumsum(var_expl)), 6),1))
if (!is.null(lm_shoot_sum)) { cat('\nShoot weight model summary (shoot_wt_eff ~ PC1 + PC2):\n'); print(lm_shoot_sum$coef); cat('\nR-squared: ', round(lm_shoot_sum$r.squared,3), '\n') }
if (!is.null(lm_florets_sum)) { cat('\nFlorets model summary (florets_eff ~ PC1 + PC2):\n'); print(lm_florets_sum$coef); cat('\nR-squared: ', round(lm_florets_sum$r.squared,3), '\n') }
cat('\nReport generated: results/report.md\n')
sink()

# Build an auto-generated markdown section (will be injected between markers)
auto_section <- c(
  '<!-- AUTO-REPORT-START -->',
  '',
  '# Auto-generated: Figures and numeric summaries',
  '',
  paste0('Generated: ', Sys.Date()),
  '',
  '## Figures (auto-inserted)',
  '',
  paste0('### 1) Drought response heatmap'),
  paste0('![](', file.path('..', plots_dir, 'drought_response_heatmap_simple.png'), ')'),
  '',
  'Interpretation: This heatmap shows standardized differences (Drought - Watered) per genotype and trait; red = increase, blue = decrease. Stars mark nominal t-test significance.',
  '',
  paste0('### 2) PCA scree plot'),
  paste0('![](', file.path('..', plots_dir, 'pca_scree_plot.png'), ')'),
  '',
  'Interpretation: The scree plot shows variance explained by each PC; PC1/PC2 often capture the major response syndromes.',
  '',
  paste0('### 3) PCA trait loadings (biplot)'),
  paste0('![](', file.path('..', plots_dir, 'pca_trait_loadings.png'), ')'),
  '',
  'Interpretation: Arrows indicate trait contributions to PC axes; similar directions = correlated drought responses.',
  '',
  paste0('### 4) PCA trait loadings colored by sensitivity'),
  paste0('![](', file.path('..', plots_dir, 'pca_trait_loadings_colored_by_sensitivity.png'), ')'),
  '',
  'Interpretation: Genotypes colored by mean effect (sensitivity) to highlight clusters.',
  '',
  paste0('### 5) PCA genotype:treatment by treatment'),
  paste0('![](', file.path('..', plots_dir, 'pca_genotype_treatment_by_treatment.png'), ')'),
  '',
  'Interpretation: Shows separation (or overlap) between drought and watered states across genotypes.',
  '',
  paste0('### 6) PCA genotype:treatment biplot'),
  paste0('![](', file.path('..', plots_dir, 'pca_genotype_treatment_biplot.png'), ')'),
  '',
  'Interpretation: Overlays trait loadings on genotype:treatment points to identify drivers of separation.',
  '',
  paste0('### 7) PCA vs Shoot weight effect'),
  paste0('![](', file.path('..', plots_dir, 'pca_shoot_weight_effect.png'), ')'),
  '',
  'Interpretation: Relationship between PC1 and genotype-level shoot weight effect; regression summary included below.',
  '',
  paste0('### 8) PCA vs Florets effect'),
  paste0('![](', file.path('..', plots_dir, 'pca_florets_effect.png'), ')'),
  '',
  'Interpretation: Relationship between PC1 and floret effect; regression summary included below.',
  '',
  '## Numeric summaries (auto-inserted)',
  '',
  paste0('- Number of significant genotype x trait t-tests (p < 0.05): **', n_significant, '**'),
  '',
  '### Top 10 absolute effects (genotype, variable, effect)',
  '',
  paste0('```', collapse=''),
  paste(capture.output(print(top_effects)), collapse='\n'),
  paste0('\n', '```'),
  '',
  '### PCA variance explained (first 6 PCs)',
  '',
  paste0('```', collapse=''),
  paste(capture.output(print(round(head(data.frame(PC = seq_along(var_expl), VarPct = 100*var_expl, CumPct = 100*cumsum(var_expl)), 6),1))), collapse='\n'),
  paste0('\n', '```')
)

# Append model summaries if present to the auto section
if (!is.null(lm_shoot_sum)) {
  auto_section <- c(auto_section, '', '### Shoot weight model: shoot_wt_eff ~ PC1 + PC2', '', '```', capture.output(print(lm_shoot_sum)), '```')
}
if (!is.null(lm_florets_sum)) {
  auto_section <- c(auto_section, '', '### Florets model: florets_eff ~ PC1 + PC2', '', '```', capture.output(print(lm_florets_sum)), '```')
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
  myco_sub <- as.data.frame(myco_sub)
  rownames(myco_sub) <- myco_sub$genotype
  myco_sub$genotype <- NULL
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
    write.csv(cor_df, file.path(results_dir, 'pca_mycorrhizal_correlations.csv'), row.names = FALSE)

    # Create scatter plots for each mycorrhizal var vs PC1
    for (var in cor_df$variable) {
      plot_df <- data.frame(genotype = rownames(myco_sub), effect = as.numeric(myco_sub[[var]]), PC1 = pc1_scores)
      p_plot <- ggplot(plot_df, aes(x = PC1, y = effect, label = genotype)) +
        geom_point() + geom_smooth(method = 'lm', se = TRUE) + ggrepel::geom_text_repel(size = 3) +
        labs(title = paste('PC1 vs', var, 'effect'), x = 'PC1', y = paste(var, 'effect (Drought - Watered)')) +
        theme_minimal()
      ggsave(file.path(plots_dir, paste0('pca_myco_', var, '_vs_PC1.png')), plot = p_plot, width = 8, height = 6, dpi = 300)
    }

    # FDR-correct p-values and pick top associations
    cor_df$padj <- p.adjust(cor_df$p.value, method = 'BH')
    # order by absolute correlation
    cor_df_ord <- cor_df %>% arrange(desc(abs(cor)))
    write.csv(cor_df_ord, file.path(results_dir, 'pca_mycorrhizal_correlations_with_fdr.csv'), row.names = FALSE)

    # top 3 associations by absolute correlation
    top_myco <- head(cor_df_ord, 3)

    # build short interpretations for the top associations
    top_lines <- apply(top_myco, 1, function(r) {
      var <- r['variable']; corv <- as.numeric(r['cor']); padj <- as.numeric(r['padj'])
      sign_txt <- ifelse(corv > 0, 'higher PC1 scores are associated with larger drought-induced increases', 'higher PC1 scores are associated with larger drought-induced decreases')
      sprintf('- %s: r=%.2f, FDR p=%.3g — %s in %s', var, corv, padj, sign_txt, var)
    })

    # add correlation table and top-variable interpretations to auto_section
    auto_section <- c(auto_section, '', '## Mycorrhizal × PC1 correlations', '', paste0('Full table written to `results/pca_mycorrhizal_correlations_with_fdr.csv`'), '')
    auto_section <- c(auto_section, '', '### Top mycorrhizal associations with PC1', '', top_lines, '')
    # also include small embedded previews (links) for the top variables
    for (v in top_myco$variable) {
      auto_section <- c(auto_section, paste0('![](', file.path('..', plots_dir, paste0('pca_myco_', v, '_vs_PC1.png')), ')'), '')
    }
  }
}

auto_section <- c(auto_section, '', '<!-- AUTO-REPORT-END -->')

# Inject auto_section into existing report.md without overwriting manual edits
start_marker <- '<!-- AUTO-REPORT-START -->'
end_marker <- '<!-- AUTO-REPORT-END -->'

if (file.exists(report_md)) {
  orig <- readLines(report_md)
  if (any(grepl(start_marker, orig)) && any(grepl(end_marker, orig))) {
    s <- which(grepl(start_marker, orig))[1]
    e <- which(grepl(end_marker, orig))[1]
    # Keep content before start_marker and after end_marker
    new_content <- c(orig[1:s], auto_section, orig[e:length(orig)])
    writeLines(new_content, con = report_md)
  } else {
    # Append auto section at the end
    new_content <- c(orig, '', auto_section)
    writeLines(new_content, con = report_md)
  }
} else {
  # Write a new report with the auto section only
  writeLines(c('# MAC-2022-Drought: Report', '', auto_section), con = report_md)
}

cat('Auto-generated section written to', report_md, '\n')
