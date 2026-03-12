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
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + ggrepel::geom_text_repel(size = 3) +
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

# ---- PC1 backing metrics ----
pc1_var_pct <- 100 * var_expl[1]

top_loadings <- data.frame(variable = rownames(pca$rotation), loading = pca$rotation[,"PC1"], stringsAsFactors = FALSE)
top_loadings <- top_loadings[order(-abs(top_loadings$loading)), ]
write.csv(top_loadings, file.path(results_dir, "pca_PC1_top_loadings_regen.csv"), row.names = FALSE)

scores_perf$PC1_group <- dplyr::ntile(scores_perf$PC1, 2)
group_summary <- scores_perf %>% dplyr::group_by(PC1_group) %>% dplyr::summarise(n = dplyr::n(), mean = mean(PC1, na.rm=TRUE), sd = sd(PC1, na.rm=TRUE), se = sd/sqrt(n), lo = mean - 1.96*se, hi = mean + 1.96*se)
write.csv(group_summary, file.path(results_dir, "pca_PC1_group_summary_regen.csv"), row.names = FALSE)

lm_group <- lm(PC1 ~ factor(PC1_group), data = scores_perf)
anova_group <- anova(lm_group)
group_r2 <- summary(lm_group)$r.squared
write.csv(data.frame(stat = c('R2'), value = c(group_r2)), file.path(results_dir, "pca_PC1_group_r2_regen.csv"), row.names = FALSE)

cohens_d <- NA_real_
if (all(c(1,2) %in% scores_perf$PC1_group)) {
  g1 <- scores_perf$PC1[scores_perf$PC1_group == 1]
  g2 <- scores_perf$PC1[scores_perf$PC1_group == 2]
  pooled_sd <- sqrt(((length(g1)-1)*var(g1, na.rm=TRUE) + (length(g2)-1)*var(g2, na.rm=TRUE)) / (length(g1)+length(g2)-2))
  cohens_d <- (mean(g2, na.rm=TRUE) - mean(g1, na.rm=TRUE)) / pooled_sd
}
write.csv(data.frame(cohens_d = cohens_d), file.path(results_dir, "pca_PC1_cohens_d_top_bottom_regen.csv"), row.names = FALSE)

perm_p <- NA_real_
if (all(c(1,2) %in% scores_perf$PC1_group)) {
  obs_diff <- mean(scores_perf$PC1[scores_perf$PC1_group==2], na.rm=TRUE) - mean(scores_perf$PC1[scores_perf$PC1_group==1], na.rm=TRUE)
  perm_diffs <- replicate(2000, {
    perm_group <- sample(scores_perf$PC1_group)
    mean(scores_perf$PC1[perm_group==2], na.rm=TRUE) - mean(scores_perf$PC1[perm_group==1], na.rm=TRUE)
  })
  perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))
}
write.csv(data.frame(perm_p = perm_p), file.path(results_dir, "pca_PC1_perm_p_top_bottom_regen.csv"), row.names = FALSE)

boot_ci <- function(x, B = 2000) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(c(NA, NA))
  sims <- replicate(B, mean(sample(x, replace = TRUE)))
  quantile(sims, c(0.025, 0.975))
}
ci1 <- boot_ci(scores_perf$PC1[scores_perf$PC1_group==1])
ci3 <- boot_ci(scores_perf$PC1[scores_perf$PC1_group==3])
write.csv(data.frame(group = c(1,3), lo = c(ci1[1], ci3[1]), hi = c(ci1[2], ci3[2])), file.path(results_dir, "pca_PC1_bootstrap_CIs_top_bottom_regen.csv"), row.names = FALSE)

if (requireNamespace('vegan', quietly = TRUE)) {
  groups_for_adonis <- scores_perf$PC1_group[match(rownames(effect_imp), scores_perf$genotype)]
  if (length(groups_for_adonis) == nrow(effect_imp)) {
    ad <- vegan::adonis2(as.matrix(effect_imp) ~ groups_for_adonis, permutations = 999)
    capture_out <- capture.output(print(ad))
    writeLines(capture_out, con = file.path(results_dir, "pca_PC1_permanova_regen.txt"))
  }
}

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
write.csv(corr_df, file.path(results_dir, 'pca_PC1_correlations_regen.csv'), row.names = FALSE)

# Append to summary txt
sink(summary_txt, append = TRUE)
cat('\n---- PC1 backing metrics (regenerated) ----\n')
cat(sprintf('PC1 variance (%%): %.2f\n', pc1_var_pct))
cat('Top loadings (top 6):\n'); print(head(top_loadings, 6))
cat('\nGroup summary (PC1 tertiles):\n'); print(group_summary)
cat('\nANOVA (PC1 ~ tertile):\n'); print(anova_group)
cat('\nR-squared (PC1 ~ tertile):\n'); print(group_r2)
cat('\nCohen\'s d (top vs bottom tertiles):\n'); print(cohens_d)
cat('\nPermutation p (top vs bottom tertiles):\n'); print(perm_p)
cat('\nCorrelation summary (PC1 vs variables):\n'); print(corr_df)
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

# ---- Trait grouping by PC1 groups (regenerated; automatic k selection) ----
set.seed(42)
# restrict maximum tested k to 3 by default
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
    try({
      sink(file.path(results_dir, 'pca_silhouette_log_regen.txt'), append = TRUE)
      cat('silhouette log for possible ks (regen):\n')
      print(sil_log)
      sink()
    }, silent = TRUE)
    if (all(is.na(sil_scores))) chosen_k <- 2 else chosen_k <- possible_k[which.max(sil_scores)]
  } else {
    chosen_k <- 2
  }
} else chosen_k <- 1
write.csv(data.frame(chosen_k = chosen_k), file.path(results_dir, 'pca_chosen_k_regen.csv'), row.names = FALSE)
trait_group_k <- chosen_k
if (nrow(scores) >= trait_group_k) {
  if (trait_group_k == 1) {
    raw_clusters <- rep(1, nrow(scores))
    cm <- mean(scores$PC1, na.rm = TRUE)
    clusters_ord <- raw_clusters
    cluster_for_geno <- clusters_ord[match(rownames(effect_imp), scores$genotype)]
  # Diagnostic logging
  dbg_file <- file.path(results_dir, paste0('pca_trait_group_debug_k', trait_group_k, '_regen.txt'))
  dir.create(dirname(dbg_file), showWarnings = FALSE, recursive = TRUE)
  sink(dbg_file)
  cat('Diagnostic for trait grouping by PC1 (regen)\n')
  cat('trait_group_k =', trait_group_k, '\n')
  cat('nrow(effect_imp) =', nrow(effect_imp), '\n')
  cat('length(scores$genotype) =', length(scores$genotype), '\n')
  cat('first 10 effect_imp rownames:\n'); print(head(rownames(effect_imp), 10))
  cat('first 10 scores genotypes:\n'); print(head(scores$genotype, 10))
  matched <- rownames(effect_imp) %in% scores$genotype
  cat('matched count:', sum(matched), 'of', length(matched), '\n')
  if (any(!matched)) { cat('unmatched (first 20):\n'); print(head(rownames(effect_imp)[!matched], 20)) }
  cat('raw_clusters length:', length(raw_clusters), '\n')
  cat('clusters_ord length:', length(clusters_ord), '\n')
  cat('cluster_for_geno length:', length(cluster_for_geno), '\n')
  cat('cluster_for_geno summary:\n'); print(summary(cluster_for_geno))
  sink()
  } else {
    raw_clusters <- kmeans(scores$PC1, centers = trait_group_k)$cluster
  cm <- tapply(scores$PC1, raw_clusters, mean)
  ord <- order(cm)
  remap <- setNames(seq_along(ord), names(cm)[ord])
  clusters_ord <- remap[as.character(raw_clusters)]
    cluster_for_geno <- clusters_ord[match(rownames(effect_imp), scores$genotype)]
  }
  trait_stats_list <- lapply(colnames(effect_imp), function(tr) {
    vals <- effect_imp[[tr]]
    # compute group means safely
    ok_idx <- !is.na(cluster_for_geno)
    group_means_full <- rep(NA_real_, trait_group_k)
    names(group_means_full) <- seq_len(trait_group_k)
    if (any(ok_idx)) {
      mm <- tapply(vals[ok_idx], cluster_for_geno[ok_idx], mean, na.rm = TRUE)
      if (length(mm) > 0) {
        group_means_full[names(mm)] <- as.numeric(mm)
      }
    }
    group_ns <- sapply(seq_len(trait_group_k), function(g) sum(!is.na(vals) & cluster_for_geno == g))
    pval <- NA_real_
    if (trait_group_k == 2) {
      if (all(group_ns >= 2, na.rm = TRUE)) {
        g1 <- vals[cluster_for_geno == 1]
        g2 <- vals[cluster_for_geno == 2]
        if (sum(!is.na(g1)) >= 2 && sum(!is.na(g2)) >= 2) {
          pval <- tryCatch(t.test(g1, g2)$p.value, error = function(e) NA_real_)
        }
      }
    } else if (trait_group_k > 2) {
      df_tv <- data.frame(val = vals, grp = cluster_for_geno)
      df_tv <- df_tv[!is.na(df_tv$val) & !is.na(df_tv$grp), , drop = FALSE]
      if (nrow(df_tv) >= 3) {
        gs <- table(df_tv$grp)
        if (sum(gs >= 2) >= 2) {
          a_mod <- tryCatch(aov(val ~ factor(grp), data = df_tv), error = function(e) NULL)
          if (!is.null(a_mod)) {
            an <- tryCatch(anova(a_mod), error = function(e) NULL)
            if (!is.null(an) && nrow(an) >= 1) {
              pval <- as.numeric(an$"Pr(>F)"[1])
            }
          }
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
  write.csv(trait_stats, file.path(results_dir, paste0('pca_trait_group_stats_PC1_k', trait_group_k, '_regen.csv')), row.names = FALSE)
  means_mat <- t(sapply(colnames(effect_imp), function(tr) {
    mm <- tapply(effect_imp[[tr]], cluster_for_geno, mean, na.rm = TRUE)
    out <- rep(NA_real_, trait_group_k)
    if (length(mm) > 0) {
      idx <- as.integer(names(mm))
      out[idx] <- as.numeric(mm)
    }
    names(out) <- paste0('group', seq_len(trait_group_k))
    out
  }))
  means_df <- data.frame(trait = rownames(means_mat), means_mat, stringsAsFactors = FALSE, row.names = NULL)
  write.csv(means_df, file.path(results_dir, paste0('pca_trait_group_means_PC1_k', trait_group_k, '_regen.csv')), row.names = FALSE)
  long_means <- tidyr::pivot_longer(means_df, cols = starts_with('group'), names_to = 'group', values_to = 'mean')
  p_heat <- ggplot(long_means, aes(x = group, y = trait, fill = mean)) + geom_tile() + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', na.value = 'grey50') + labs(title = paste('Trait means by PC1 groups (k=', trait_group_k, ')', sep='')) + theme_minimal()
  ggsave(file.path(plots_dir, paste0('pca_PC1_trait_group_means_k', trait_group_k, '_regen.png')), plot = p_heat, width = 8, height = max(6, 0.25 * nrow(means_df)), dpi = 300)
  sink(summary_txt, append = TRUE)
  cat('\n---- Trait grouping by PC1 (k=', trait_group_k, ') (regenerated) ----\n', sep='')
  cat('Trait group counts:\n')
  print(table(trait_stats$assigned_group))
  cat('\nTop traits assigned to each group (by mean):\n')
  for (g in seq_len(trait_group_k)) {
    topg <- means_df %>% dplyr::arrange(dplyr::desc(.data[[paste0('group', g)]])) %>% dplyr::slice_head(n = 8) %>% dplyr::pull(trait)
    cat('Group', g, ':', paste(topg, collapse = ', '), '\n')
  }
  sink()
}

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
        geom_point() + geom_smooth(method = 'lm', formula = y ~ x, se = TRUE) + ggrepel::geom_text_repel(size = 3) +
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
