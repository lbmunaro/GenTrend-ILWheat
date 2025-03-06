# R script used in Munaro et al. 2025
#
# Script author: Lucas Berger Munaro

# Script objectives ----
# - Organize and clean phenotypic data from advanced yield trials
# - Run BLUP models to estimate reliability of each trait within each trial
# - Run BLUE models to estimate genotype's blues and standard error for each trait within each trial

rm(list=objects()) # Clear workspace

# Packages ----
library(tidyverse) # R packages for data science
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data
library(asreml) # ASReml-R package

# Phenotypic data ----

# Read, clean, and organize data from the Breedbase
pheno_w <- # Phenotypic data in wide format
  read.csv('Data/GenTrend-ILWheat_data_2022nov10.csv', skip = 3) |> 
  remove_empty(which = c('cols')) |> # Remove any completely empty columns
  clean_names() |> # Standardize column names for consistency
  # Select variables relevant to the analysis
  dplyr::select(observation_unit_name,study_name,study_year,location_name,
                block_number, row_number, col_number, germplasm_name,
                heading_time_julian_date_jd_co_321_0001233,
                plant_height_cm_co_321_0001301,
                grain_yield_kg_ha_co_321_0001218,
                grain_test_weight_g_l_co_321_0001210) |>
  # Rename variables for easier reference
  rename(id=observation_unit_name,
         year=study_year,
         loc=location_name,
         env=study_name,
         blk=block_number,
         row=row_number,
         col=col_number,
         gen=germplasm_name,
         heading_time=heading_time_julian_date_jd_co_321_0001233,
         plant_height=plant_height_cm_co_321_0001301,
         grain_yield=grain_yield_kg_ha_co_321_0001218,
         test_weight=grain_test_weight_g_l_co_321_0001210) |>
  arrange(year,loc,gen,blk,row,col) |> # Arrange data by trial structure
  # Standardize location names
  mutate(loc = str_replace(loc,'St. Jacob Township, IL','St. Jacob, IL')) |>
  # Convert to factors
  mutate_at(vars(id:gen),as.factor) |>
  # Convert to numeric format
  mutate_at(vars(heading_time:test_weight),as.numeric) |>
  # Replace zero values with NA
  mutate_at(vars(heading_time:test_weight), ~ifelse(.==0,NA,.)) |>
  # Exclude Advanced High Yield (AdvHY) trials
  filter(!grepl('AdvHY',env)) |>
  # Drop unused levels
  droplevels() |>
  glimpse()

# Convert to long format
pheno <- pheno_w |>
  pivot_longer(cols = heading_time:test_weight,
               names_to = 'var', # Trait names
               values_to = 'val', # Trait values
               values_drop_na = T) |>
  mutate(var=as.factor(var)) |>
  glimpse()

# Single-trial analysis ----

## BLUP model ----
# Function to calculate BLUPs and reliability of each trait within each trial
mod.blups <- function(dat){
  yr <- unique(as.character(dat$year)) # Identify the year being analyzed
  if('2021' %in% yr){
    # Model structure differs for 2021, which includes spatial row/column adjustment
    mod <- asreml(fixed = val~1+blk,
                  random = ~gen+row+col,
                  data = dat, na.action = na.method(y='include', x='include'))
  }else{
    mod <- asreml(fixed = val~1+blk,
                  random = ~gen,
                  data = dat, na.action = na.method(y='include', x='include'))
  }
  
  # Update model
  mod <- update(mod)
  
  # Extract BLUPs and calculate reliability
  blups<- predict(mod, classify='gen', ignore=c('(Intercept)','blk'))$pvals
  pev<- blups[,'std.error']^2 # Prediction error variance
  Vg<- summary(mod)$varcomp['gen','component'] # Genetic variance estimate
  rel<- 1-(pev/Vg) # Reliability for each genotype
  m_rel<- mean(rel) # Mean reliability
  
  # Return mean reliability
  res<- data.frame(m_rel=m_rel)
  return(res)
}

### Run BLUP model ----
blups <- pheno|>
  group_by(env,var) |> 
  # Filter out environments with only one block, as BLUPs require replication
  mutate(no_blk=length(unique(blk))) |> 
  filter(no_blk>1) |> 
  dplyr::select(-no_blk) |>
  nest() |>
  # Apply BLUP function for each combination of trial and trait
  mutate(blups=map(data, ~mod.blups(.x))) |>
  unnest_wider(blups) |>
  ungroup() |>
  mutate_if(is.numeric,~round(.,3)) |>
  glimpse()

### Summarize reliability results ----
blups |>
  dplyr::select(-c(data)) |>
  group_by(var) |>
  summarise(min=min(m_rel),
            max=max(m_rel),
            mean=mean(m_rel),
            median=(median(m_rel))) |>
  mutate_if(is.numeric, ~round(.,2))

## BLUE model ----
# Function to estimate genotype's blues and standard error for each trait within each trial
mod.blues <- function(dat){
  yr <- unique(as.character(dat$year)) # Identify the year being analyzed
  
  # Model structure differs for 2021 due to spatial row/column adjustment
  if('2021' %in% yr){
    mod <- asreml(fixed = val~1+gen+blk,
                  random = ~row+col,
                  data = dat,
                  na.action = na.method(y='include', x='include'))
  }else{
    mod <- asreml(fixed = val~1+gen+blk,
                  data = dat,
                  na.action = na.method(y='include', x='include'))
  }
  
  # Update model
  mod <- update(mod)
  
  # Extract BLUEs
  blues<- predict(mod, classify='gen')$pvals
  return(blues)
}

### Run BLUE model ----
out_blues <- blups |>
  unnest(data) |>
  group_by(env,var) |>
  nest() |>
  # Apply BLUE function for each combination of trial and trait
  mutate(blues=map(data, ~mod.blues(.x))) |>
  dplyr::select(-data) |>
  unnest(blues) |>
  clean_names() |>
  rename(prd_val=predicted_value) |>
  dplyr::select(-status) |>
  # Add year and location information
  left_join(pheno |>
              group_by(env) |>
              summarise(year=unique(year),
                        loc=unique(loc)),
            by='env') |>
  glimpse()

# Test weight imputation ----
# Trials with test weight data from only one replication require imputation

# Average standard error for test weight standard error from 2001 to 2013 for imputation
tw.std_error <- out_blues |>
  ungroup() |>
  filter(var=='test_weight'&as.numeric(as.character(year))<2014) |>
  summarise(se=sqrt(mean((std_error^2)*3))) |>
  unlist(use.names = F)

# Create final dataset with BLUEs and imputed values
blues <- out_blues |> 
  bind_rows(pheno|>
              group_by(env,var) |>
              mutate(no_blk=length(unique(blk))) |> 
              filter(no_blk==1) |>
              mutate(prd_val=val,
                     std_error=tw.std_error) |>
              dplyr::select(env,var,gen,prd_val,std_error,year,loc)) |>
  group_by(gen) |>
  mutate(gen=as.factor(gen),
         gidyr=min(as.numeric(as.character(year))), # First year genotype entered trial
         # First year (2001) we consider 'IL-98' lines as first entry
         gidyr=ifelse(gidyr!=2001,gidyr,ifelse(grepl('^98',gen),2001,NA)), 
         # Add NA to gidyr for genotypes we exclude from the analysis
         gidyr=ifelse(!grepl('^\\d',gen),NA,gidyr), # Do not starts with number (not IL)
         gidyr=ifelse(grepl('B-B',gen),NA,gidyr), # Not IL
         gidyr=ifelse(grepl('^96',gen),NA,gidyr), # First entered before 2001
         gidyr=ifelse(grepl('^97',gen),NA,gidyr), # First entered before 2001
         #IL genotypes that starts with US
         gidyr=ifelse(grepl('^US',gen),
                      min(as.numeric(as.character(year))),gidyr), 
         ## Gidyr for Kaskaskia = 1993, based on the line name (90-7514)
         gidyr=ifelse(gen=='Kaskaskia',1993,gidyr), 
         gidyr=as.factor(gidyr),
         year_n=as.numeric(as.character(year)),
         gidyr_n=as.numeric(as.character(gidyr)),
         wt=(1/(std_error^2)))|> # Calculate weights to include in following models
  ungroup() |>
  dplyr::select(env,year,year_n,loc,gen,gidyr,gidyr_n,var,prd_val,std_error,wt) |>
  arrange(year,loc,gen,var) |>
  droplevels.data.frame()

# Save results ----
save(pheno, blups, blues, file = 'Data/Pheno-Blups-Blues.RData')

# End ----