# R script used in Munaro et al. 2025
#
# Script author: Lucas Berger Munaro

# Script objective ----
# - Calculate the genetic correlations among traits

rm(list=objects()) # Clear workspace

# Packages ----
library(tidyverse) # R packages for data science
library(asreml) # ASReml-R package

# Data ----

## Load BLUEs ----
load('Data/pheno-blups-blues.RData')

# Prepare dataset for genetic correlation analysis
# This step fills in missing trait values with NA across the whole dataset
blues.na <- blues |>
  dplyr::select(-c(std_error, wt)) |> # Remove columns not needed for wide format
  pivot_wider(values_from = prd_val, names_from = var) |> # Pivot wider to include NA for unobserved traits
  pivot_longer(names_to = 'var', values_to = 'prd_val',
               cols = c(grain_yield:heading_time), 
               values_drop_na = FALSE) |> # Pivot longer without dropping NA values
  left_join(blues) |> # Add back standard errors and weights for analysis
  mutate_if(is.character, as.factor) |>
  arrange(var, env, gen) |>
  mutate(check = ifelse(gen == 'Kaskaskia', TRUE, FALSE), # Label check genotypes
         check_f = as.factor(check)) |> # Convert check to a factor
  filter(!is.na(gidyr) | check == TRUE) |> # Filter out genotypes without entry year (gidyr)
  droplevels() |>
  glimpse()

# Select trials where all traits (grain_yield, heading_time, plant_height, test_weight) are observed
trials <- blues |>
  # Filter to retain rows with all trait values present
  group_by(env) |>
  filter(all(c('grain_yield', 'heading_time', 'plant_height', 'test_weight') %in% unique(var)) & !is.na(prd_val)) |>
  summarise(env = unique(env)) |>
  droplevels() |>
  as.vector()

dat <- blues.na |>
  filter(env%in%trials$env) |> # Filter to include only selected trials
  dplyr::filter(!is.na(prd_val)) |>
  droplevels() |>
  arrange(var,gen,env) |>
  glimpse()

# Further organize the dataset to ensure all traits (yield, height, test weight, heading time) are non-NA
dat <- dat |>
  select(-c(std_error,wt)) |> # Wide format to check trait availability
  pivot_wider(names_from = var, values_from = prd_val) |>
  filter_at(.vars = c('grain_yield','heading_time','plant_height','test_weight'),
            ~!is.na(.)) |> # Keep only rows with all traits observed
  pivot_longer(cols = c('grain_yield', 'heading_time', 'plant_height', 'test_weight'),
               names_to = 'var', values_to = 'prd_val') |> # Return to long format
  left_join(dat) |>
  mutate(var=as.factor(var)) |>
  glimpse()

# Multi-trait model for first 5 years (2001-2005)
mt_mod_1st5yr <- asreml(fixed = prd_val~var, # Multi-trait model with each trait as a fixed effect
                        random = ~corgh(var):gen + diag(var):env + diag(var):env:gen,
                        family=asr_gaussian(dispersion = 1),
                        residual= ~units:var, 
                        na.action = na.method(y='include', x='include'), 
                        asmv=var, weights=wt, workspace='1gb',
                        data=dat|>filter(year_n <2006)) # Filter for the first 5 years
mt_mod_1st5yr <- update(mt_mod_1st5yr)

# Extract genetic correlations among traits from first 5-year model
corgh_1st5yr <- summary(mt_mod_1st5yr)$varcomp |>
  as.data.frame() |>
  rownames_to_column() |>
  mutate(component=round(component,2)) |>
  filter(str_detect(rowname, c('.cor'))) |>
  # Clean row names to retain only the trait names in the correlation
  mutate(rowname=str_remove(rowname,'var:gen!')) |>
  mutate(rowname=str_remove(rowname,'var!')) |>
  mutate(rowname=str_remove(rowname,'!var!')) |>
  mutate(rowname=str_remove(rowname,'.cor')) |>
  mutate(df='1st5yr') |>
  glimpse()

# Multi-trait model for last 5 years (2017-2021)
mt_mod_lst5yr <- asreml(fixed = prd_val~var,
                        random = ~corgh(var):gen + diag(var):env + diag(var):env:gen,
                        family=asr_gaussian(dispersion = 1),
                        residual= ~units:var, 
                        na.action = na.method(y='include', x='include'), 
                        asmv=var, weights=wt, workspace='1gb',
                        data=dat|>filter(year_n >2016)) # Filter for the last 5 years
mt_mod_lst5yr <- update(mt_mod_lst5yr)

# Extract genetic correlations among traits from last 5-year model
corgh_lst5yr <- summary(mt_mod_lst5yr)$varcomp |>
  as.data.frame() |>
  rownames_to_column() |>
  mutate(component=round(component,2)) |>
  filter(str_detect(rowname, c('.cor'))) |>
  # Clean row names for better readability
  mutate(rowname=str_remove(rowname,'var:gen!')) |>
  mutate(rowname=str_remove(rowname,'var!')) |>
  mutate(rowname=str_remove(rowname,'!var!')) |>
  mutate(rowname=str_remove(rowname,'.cor')) |>
  mutate(df='lst5yr') |>
  glimpse()

# Combine correlation data across time periods
print(corgh_1st5yr)
print(corgh_lst5yr)

gencorr <- rbind(corgh_1st5yr, corgh_lst5yr) |>
  mutate(std.error=round(std.error,2)) |>
  # Format genetic correlation with standard error
  mutate(gencorr=paste(component,std.error,sep='±')) |>
  select(rowname,gencorr,df) |>
  pivot_wider(names_from = df, values_from = gencorr) |>
  rename(correlation=rowname) |>
  glimpse()
gencorr

# Save results ----
save(gencorr, file='Data/GenCorr.Rdata')

# End ----