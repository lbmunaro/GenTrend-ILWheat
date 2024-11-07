# GenTrend-ILWheat

## Overview
This repository contains R scripts developed for analyzing genetic, non-genetic, and phenotypic trends in winter wheat, spanning 21 years of breeding data from the University of Illinois. This analysis is detailed in the manuscript:  
**Munaro et al., 2025. "Genetic Trends Due to 21 Years of Winter Wheat Breeding at the University of Illinois."**

The repository is organized to streamline the analysis, from data preparation to trend analysis and genetic correlation among traits.

## Repository Contents

### Scripts

1. **01_Data_Preparation_Single_Trial_Analysis.R**  
   - Prepares and cleans data, calculates single-trial BLUPs (Best Linear Unbiased Predictions) and BLUEs (Best Linear Unbiased Estimates) for each genotype across multiple environments.

2. **02_Gen_Nongen_Pheno_Trend_Analysis.R**  
   - Conducts trend analysis to estimate genetic, non-genetic (environmental), and total phenotypic trends for each trait over the study period.

3. **03_Genetic_Correlation_Analysis.R**  
   - Computes genetic correlations among traits (grain yield, heading time, plant height, and test weight).

4. **04_Tables_Figures.R**  
   - Generates tables and figures present in the manuscript.

### Data
- **Data/**: This directory contains data files, including
  - `GenTrend-ILWheat_data_2022nov10.csv`: Phenotypic data from the University of Illinois wheat breeding program's advanced yield trials from 2001 to 2021. The dataset is available on The Triticae Toolbox (T3) database accessible at: https://wheat.triticeaetoolbox.org/.
  - `Pheno-Blups-Blues.RData`: Dataset with genotype's blues and weights for each trait within each trial, single-trial reliabilities, and clean phenotypic data.
  - `Trends.RData`: Dataset with results from the genetic, nongenetic, and phenotypic trend analyses.
  - `GenCorr.RData`: Dataset with results from the genetic, nongenetic, and phenotypic trend analyses.
  - `2022_fsa_acres_web_082222.csv`: 2022 USDA-FSA crop acreage data reported used to create Figure 2. This data was obtained from https://www.fsa.usda.gov/tools/informational/freedom-information-act-foia/electronic-reading-room/frequently-requested/crop-acreage-data

## Usage
Each script is designed to be run in sequence, from `01_` to `04_`. The analysis requires the ASReml-R package for mixed model analysis.
