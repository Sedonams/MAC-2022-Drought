## Sorghum MAC-2022 Drought Project

**Last updated:** August 5, 2025

This repository contains scripts and data for the Sorghum MAC-2022 drought experiment. The goal is to analyze drought tolerance and related traits in sorghum genotypes using a standardized workflow.

### Project Structure

- `data/` — Raw and cleaned datasets (Excel and csv files)
- `scripts/` — R scripts for data cleaning, analysis, and visualization

### Getting Started

1. Clone or download this repository.
2. Open the R scripts in the `scripts/` folder.
3. Run the cleanup script (`MAC22_data_prep.R`) to prepare the dataset for analysis. This script standardizes column names, data types, and removes inconsistencies.
4. Work on some analyses. Sedona has hers in `MAC 22 (1).R`, Bea is working on some things based off Sedona's script in `MAC22_BB.R`. 

### Main Scripts

- `MAC22_data_prep.R`: Main cleanup script. Use this to generate a uniform dataset (located in `/data` but feel free to run again).
- Other scripts in the folder provide additional analyses and visualizations.

### Data

All data files are located in the `data/` folder. The main dataset is `MAC 2022 Mother for R.xlsx`, and the cleaned up version is `MAC22_cleaned.csv`.

### How to Contribute

If you have improvements or new analyses, please add your scripts to the `scripts/` folder and update this README with instructions.

---

## Analysis scripts (detailed)

This section documents the key analysis scripts, how to run them, what they produce, and useful options.

- `scripts/MAC22_BB.R` (main analysis)
	- Purpose: Full analysis pipeline. Builds the genotype × trait drought-effect matrix (genotype mean(Drought) - mean(Watered)), runs PCA, creates plots (scree, loadings, genotype biplots), tests associations between PC1 and performance traits, computes correlations between PC1 and mycorrhizal variables, and groups traits by PC1 means.
	- Run: `Rscript -e "source('scripts/MAC22_BB.R')"`
	- Outputs (examples):
		- `results/pca_chosen_k.csv` — chosen k for trait grouping (silhouette-based by default)
		- `results/pca_trait_group_stats_PC1_k{K}.csv` — per-trait assigned group, p.value and BH-adjusted padj
		- `results/pca_trait_group_means_PC1_k{K}.csv` — means per group
		- `results/pca_trait_group_debug_k{K}.txt` — diagnostic log (matching, cluster sizes)
		- `plots/pca_trait_loadings.png` and `plots/pca_trait_loadings_colored_by_sensitivity.png` — biplot images
		- `plots/pca_trait_loadings_with_group_ellipses.png` — biplot with ellipses for detected groups
		- `plots/pca_trait_loadings_with_2group_ellipses.png` — visual-only 2-group ellipse biplot

- `scripts/generate_report.R` (report generator)
	- Purpose: Recompute selected outputs and insert an auto-generated block into `results/report.md`. Useful for regenerating key tables/plots without rerunning the full pipeline.
	- Run: `Rscript -e "source('scripts/generate_report.R')"`

## Trait grouping behavior and notes

- Automatic k selection: the script evaluates candidate k values using silhouette scores and selects the best k by max silhouette. For small datasets, the script restricts the maximum k to avoid tiny clusters.
- Assigned group: for each trait, the script computes the mean effect per detected group and assigns the trait to the group with the largest mean (which.max). Therefore, even with k=3 you may see only groups 1 and 3 listed as `assigned_group` for many traits if group2 was never the maximum for those traits.
- Statistical tests: for k=2 the script runs a t-test (if each group has >=2 observations). For k>2 it tries ANOVA (fallback to Kruskal–Wallis if needed). p.values are set to NA if groups are too small or test fails.

## Quick fixes and options

- Force two groups: to force two-group comparisons, add `chosen_k <- 2` after the code that computes `chosen_k` in `scripts/MAC22_BB.R`.
- Merge tiny cluster: if one cluster is tiny, you can request an automatic merge (reassign its genotypes to the nearest cluster by PC1 distance). Contact me and I can add that behavior.
- Ellipses: the script draws ellipses using `stat_ellipse(..., geom='polygon')` at `level = 0.68` (approx. 1 SD). Change `level` to `0.95` for a 95% ellipse.

## Troubleshooting

- Many NA p.values: open `results/pca_trait_group_debug_k{K}.txt` to inspect cluster sizes and `cluster_for_geno` mapping. NA p.values usually mean groups had insufficient non-missing data for the test.
- Mismatched genotype names: ensure `effect_imp` rownames and PCA `scores$genotype` match; the debug file shows the matching statistics.


---
Generated: 2025-08-18

