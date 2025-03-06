# R script used in Munaro et al. 2025
#
# Script author: Lucas Berger Munaro

# Script objective ----
# - Analyze genetic, nongenetic, and phenotypic trends over 17 years with two checks

rm(list=objects()) # Clear workspace

# Packages ----
library(tidyverse) # R packages for data science
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data
library(asreml) # ASReml-R package

# Data ----

## Load BLUEs from previous analysis ----
load('Data/pheno-blups-blues.RData')

dat_trends_17yr <-
  blues |> # BLUEs from single-trial analysis
  filter(year%in%c(as.character(2005:2021))) %>%
  mutate(check=ifelse(gen%in%c("Kaskaskia","02-18228"),T,F), # Label check genotypes
         check_f=as.factor(check)) |> # Convert check label to a factor
  filter(!is.na(gidyr)|check==T) |> # Keep only genotypes with known entry years (gidyr)
  filter(check==T|year_n==gidyr_n) |> # For entry genotypes, retain only data from the first test year
  glimpse()

# Models ----

# Function to run the genetic, nongenetic, and phenotypic trends
mod.trends_17yr <- 
  function(dat){
    var <- unique(as.character(dat$variable)) # Variable being analyzed
    if('heading_time' %in% var) { # Specific model structure for heading time
      
      # Genetic trend: Estimates genetic trend over entry year (gidyr_n) with checks
      mod_gen <- asreml(fixed=prd_val ~gidyr_n + check_f,
                        random=~year + at(year,'2021'):gen + at(year,'2021'):loc,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'), data=dat)
      mod_gen <- update(mod_gen)
      
      # Fixed effects for visualization of genetic trend (as point estimates)
      mod_gen_pt <- asreml(fixed=prd_val ~gidyr + check_f,
                           random=~year+at(year,'2021'):loc,
                           weights=wt, family=asr_gaussian(dispersion=1),
                           na.action=na.method('omit'), data=dat)
      mod_gen_pt <- update(mod_gen_pt)
      
      # Nongenetic (environmental) trend: Only check genotypes over years
      mod_env <- asreml(fixed=prd_val ~year_n,
                        random=~at(year,'2021'):loc + gen,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'),data=dat|>filter(check==T))
      mod_env <- update(mod_env)
      
      # Phenotypic trend: Only entry genotypes over years
      mod_pheno <- asreml(fixed=prd_val ~gidyr_n,
                          random=~at(year,'2021'):gen+at(year,'2021'):loc,
                          weights=wt, family=asr_gaussian(dispersion=1),
                          na.action=na.method('omit'),data=dat|>filter(check==F))
      mod_pheno <- update(mod_pheno)
      
    }else{ # Model structure for grain yield, test weight, and plant height
      
      # Genetic trend for other traits
      mod_gen <- asreml(fixed=prd_val~ gidyr_n+ check_f + loc + check_f:loc,
                        random=~ year + gen:year + gen:loc,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'), data=dat)
      mod_gen <- update(mod_gen)
      
      # Fixed effects for visualization of genetic trend (as point estimates)
      mod_gen_pt <- asreml(fixed=prd_val~ gidyr+ check_f + loc + check_f:loc,
                           random=~ year + gen:year + gen:loc,
                           weights=wt, family=asr_gaussian(dispersion=1),
                           na.action=na.method('omit'), data=dat)
      mod_gen_pt <- update(mod_gen_pt)
      
      # Nongenetic (environmental) trend with checks over years
      mod_env <- asreml(fixed=prd_val~ year_n, 
                        random=~ gen:loc:year + loc + year,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'),
                        data=dat|>filter(check==T))
      mod_env <- update(mod_env)
      
      # Phenotypic trend for entry genotypes only
      mod_pheno <- asreml(fixed=prd_val~ gidyr_n + loc, 
                          random=~ year + gen:loc + gen:year,
                          weights=wt, family=asr_gaussian(dispersion=1),
                          na.action=na.method('omit'),
                          data=dat|>filter(check==F))
      mod_pheno <- update(mod_pheno)
    }
    
    # Genetic trend figure data
    ## Predict genetic trend line for visualization
    gen_line <- predict(mod_gen,  classify='gidyr_n:check_f',
                        levels=list(check_f="FALSE", gidyr_n=c(unique(dat$year_n))))$pvals |>
      as.data.frame() |> clean_names() |> 
      rename(pval_line=predicted_value, se_line=std_error) |>
      mutate(gidyr=factor(gidyr_n),
             ci1=pval_line+1.96*se_line, ci2=pval_line-1.96*se_line)
    
    ## Predict genetic trend points for specific years
    gen_pts <- predict(mod_gen_pt, classify='gidyr:check_f',
                       levels=list(check_f="FALSE", gidyr=factor(c(unique(dat$year_n)))))$pvals |>
      as.data.frame() |>clean_names() |> 
      rename(pval_point=predicted_value,se_point=std_error)
    
    ## Combine trend line and point predictions
    dat_gen_trend <- left_join(gen_line, gen_pts) |> nest(dat_gen_trend = everything())
    
    # Wald test of fixed effects to assess significance
    wald_df <- function(mod,trend){
      aov_df <- wald(mod) |> as.data.frame() |> rownames_to_column('term') |> 
        mutate(trend=trend) |> relocate(trend)
      return(aov_df)
    }
    wald <- bind_rows(wald_df(mod_gen,'gen'), 
                      wald_df(mod_env,'env'), 
                      wald_df(mod_pheno,'pheno'))
    
    # Fixed effects estimates
    fixeff_df <- function(mod,trend){
      fixeff_df <- summary(mod, coef=TRUE)$coef.fixed |> as.data.frame() |> 
        rownames_to_column('term') |> mutate(trend=trend) |> relocate(trend)
      return(fixeff_df)
    }
    coef_fixeff <- bind_rows(fixeff_df(mod_gen,'gen'), 
                             fixeff_df(mod_env,'env'), 
                             fixeff_df(mod_pheno,'pheno'))
    raw_fixeff <- bind_rows(list(wald=wald,coef_fixeff=coef_fixeff),.id='id') |>
      clean_names()|>nest(raw_fixeff=everything())
    
    ## Summarized fixed effects
    fixeff <- full_join(wald,coef_fixeff,by=c('trend', 'term')) |> 
      clean_names() |> filter(!grepl('Intercept|residual',term)) |> 
      dplyr::select(-c(sum_of_sq,z_ratio)) |> nest(fixeff=everything())
    
    # Variance components
    varcomp_df <- function(mod,trend){
      varcomp_df <- summary(mod)$varcomp |> as.data.frame() |> 
        rownames_to_column('term') |> mutate(trend=trend) |> relocate(trend)
      return(varcomp_df)
    }
    varcomp <- bind_rows(varcomp_df(mod_gen,'gen'), 
                         varcomp_df(mod_env,'env'), 
                         varcomp_df(mod_pheno,'pheno')) |>
      clean_names() |> nest(varcomp=everything())
    
    return(data.frame(dat_gen_trend=dat_gen_trend,
                      fixeff=fixeff,
                      varcomp=varcomp))
  }

## Run model for each trait ----
out_trends_17yr <- dat_trends_17yr |>
  mutate(variable=as.character(var)) |>
  group_by(var) |>
  nest() |>
  mutate(trends=map(data,
                    ~mod.trends_17yr(.x))) |>
  dplyr::select(-data) |>
  unnest(trends) |>
  ungroup()

out_trends_17yr |> glimpse()

# Save results ----
save(out_trends_17yr, file = 'data/Trends_17yr.RData')

# End ----