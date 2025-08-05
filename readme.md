## Sorghum MAC-2022 Drought Project

**Last updated:** August 5, 2025

This repository contains scripts and data for the Sorghum MAC-2022 drought experiment. The goal is to analyze drought tolerance and related traits in sorghum genotypes using a standardized workflow.

### Project Structure

- `data/` — Raw and cleaned datasets (Excel files)
- `scripts/` — R scripts for data cleaning, analysis, and visualization

### Getting Started

1. Clone or download this repository.
2. Open the R scripts in the `scripts/` folder.
3. Run the cleanup script (`MAC22_data_prep.R`) to prepare the dataset for analysis. This script standardizes column names, data types, and removes inconsistencies.
4. Work on some analyses. Sedona has hers in `MAC 22 (1).R`, Bea is working on some things based off Sedona's script in `MAC22_BB.R`. 

### Main Scripts

- `MAC22_BB.R`: Main cleanup script. Use this to generate a uniform dataset (located in `/data` but feel free to run again).
- Other scripts in the folder provide additional analyses and visualizations.

### Data

All data files are located in the `data/` folder. The main dataset is `MAC 2022 Mother for R.xlsx`, and the cleaned up version is `MAC22_cleaned.csv`.

### How to Contribute

If you have improvements or new analyses, please add your scripts to the `scripts/` folder and update this README with instructions.

