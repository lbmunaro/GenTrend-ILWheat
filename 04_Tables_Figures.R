# R script used in Munaro et al. 2025
#
# Script author: Lucas Berger Munaro

# Script objective ----
# - Create tables and figures

rm(list=objects()) # Clear workspace

# Packages ----
library(tidyverse) # R packages for data science
library(ggspatial) # Spatial Data Framework for ggplot2
library(RColorBrewer) # ColorBrewer Palettes
library(ggtext) # Improved text rendering support for ggplot2
library(patchwork) # The Composer of Plots

# Load results ----
load('Data/Pheno-Blups-Blues.RData')
load('Data/Trends.RData')
load('Data/GenCorr.Rdata')
load('Data/Trends_17yr.RData')

# Tables ----

## Table 1 ----
# Rates of change in genetic, nongenetic, and total phenotypic values for soft red 
# winter wheat from the University of Illinois wheat breeding program from 2001 to 2021 for
# grain yield, plant height, test weight, and heading time

out_trends |>
  dplyr::select(var, fixeff) |>
  unnest(fixeff) |>
  filter(term%in%c('year_n','gidyr_n','at(check_f, FALSE):gidyr_n')) |>
  mutate(slope=round(solution,2),
         se=round(std_error,2),
         sign=ifelse(pr_chisq<0.001,'***',
                     ifelse(pr_chisq<0.01,'**',
                            ifelse(pr_chisq<0.05,'*','')))) |>
  mutate(value=paste(slope, '±', se, sign)) |>
  dplyr::select(var,trend,value) |>
  pivot_wider(names_from = trend, values_from = value)

## Table 2 ----
# Genetic correlation among traits across the first five years (2001 to 2005) and the
# last five years (2017 to 2021). The standard error follows the genetic correlation 
# coefficients
gencorr

# Figures ----

## Figure 2 ----
# Illinois map showing the six test locations. The orange shades represent the 2022 
# wheat planted hectares per county (USDA-FSA, 2022)

loc_coord <- data.frame(loc=c('Brownstown', 'Carmi', 'Neoga', 'Ridgway', 'St. Jacob', 'Urbana'),
                        lat=c(38.9948827,38.08553,39.23274,37.79930,38.74573,40.05833),
                        long=c(-88.95961,-88.19441,-88.38207,-88.27799,-89.78045,-88.22937)) |>
  glimpse()

states <- c('illinois')
highlight_states <- c('illinois')

map_data <- map_data('state') %>%
  filter(region %in% states)

county <- map_data('county') |>
  filter(region %in% states) |>
  mutate(subregion=str_replace_all(subregion,' ',''),
         subregion=str_replace_all(subregion,'[.]',''))

wheat_ha <- read.csv('Data/2022_fsa_acres_web_082222.csv') |>
  mutate_if(is.character,~tolower(.)) |>
  dplyr::filter(Crop%in%'wheat') |>
  mutate(Planted.Acres=str_replace(Planted.Acres,',',''),
         Planted.Acres=as.numeric(Planted.Acres)) |>
  group_by(State,County,Crop) |>
  summarise(Planted.Acres=sum(Planted.Acres,na.rm=T)) |>
  mutate(region=State,
         subregion=County,
         subregion=str_replace_all(subregion,
                                   c('dewitt'='de witt','dupage'='du page','st. clair'='st clair','dekalb'='de kalb' ))) |>
  mutate(subregion=str_replace_all(subregion,' ',''),
         subregion=str_replace_all(subregion,'[.]',''))|>
  mutate(Planted.Hectares = round(Planted.Acres*0.404686,0)) |>
  glimpse()

breaks <- c(1, 500, 1000, 2500, 5000, 10000, 20000, 40000)
labels <- c('< 500', '501 to 1,000', '1,001 to 2,500', '2,501 to 5,000',
            '5,001 to 10,000', '10,001 to 20,000', '> 20,000')

county_ha <- county |>
  left_join(wheat_ha) |>
  mutate(Planted.Hectares = ifelse(Planted.Hectares == 0, NA, Planted.Hectares)) |>
  mutate(classes = cut(Planted.Hectares, breaks=breaks, labels = labels, include.lowest = TRUE, right = FALSE)) |>
  glimpse()

ggplot() +
  geom_polygon(data = county_ha, aes(x = long, y = lat, group = group, fill = classes), color = 'gray70') +
  scale_fill_brewer(name = 'Planted area (ha)', palette = 'Oranges', na.value = 'white', drop= TRUE, 
                    labels = c('< 500', '501 to 1,000', '1,001 to 2,500', '2,501 to 5,000',
                               '5,001 to 10,000', '10,001 to 20,000', '> 20,000',
                               'Not estimated')) +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), fill = 'transparent', color = 'gray50') +
  geom_point(data = loc_coord, aes(x = long, y = lat, group = NA), color = '#13294B', size = 1) +
  geom_point(data = loc_coord, aes(x = long, y = lat), color = '#13294B', size = 2, shape = 21) + 
  geom_text(data = loc_coord, aes(x = long, y = lat, group = NA, label = loc),
            angle = 0, vjust = 0.5, hjust = -0.1, color = '#13294B', size = 4) +
  coord_fixed(1.3) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.background=element_rect(fill='white'),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.key=element_rect(fill='white',color=NA),
        legend.position = 'left',
        plot.margin = margin(0,0,0,0, 'in')) +
  xlim(range(loc_coord$long) + c(-2,1)) +
  annotation_north_arrow(location='bl',pad_x=unit(0,'in'),pad_y=unit(0.2,'in'),
                         style = north_arrow_nautical,
                         height = unit(0.5,'in'),width = unit(0.5,'in'))
ggsave('Figures/Figure2.tif', width = 7, height = 5, units = 'in', dpi = 300)

## Figure 3 ----
# Traits evaluated within each trial (combination of location and year). The colors 
# represent the number of genotypes evaluated within each trial, and the white tile represents 
# an unobserved trait.

fig3_summary <- pheno |>
  group_by(var) |>
  summarise(nyear = length(unique(year)),
            nloc = length(unique(loc)),
            nenv = length(unique(env)),
            ngen = length(unique(gen)),
            tot_obs = n()) |>
  mutate(tot_obs = formatC(tot_obs, format = 'f', big.mark = ',', digits = 0),
         nenv = formatC(nenv, format = 'f', big.mark = ',', digits = 0))

pheno |>
  group_by(year, loc, var) |>
  summarise(no_gen = length(unique(gen)), .groups = 'drop') |>
  pivot_wider(names_from = loc, values_from = no_gen) |>
  pivot_longer(names_to = 'loc', values_to = 'no_gen', cols = 3:8) |>
  ggplot(aes(year, loc)) +
  geom_tile(aes(fill = no_gen), color = 'black', na.rm = TRUE) +
  facet_wrap(~var, labeller = labeller(var = c(
    'grain_yield' = paste0('Grain yield\n(Trials: ', fig3_summary$nenv[fig3_summary$var == 'grain_yield'], 
                           ', Obs.: ', fig3_summary$tot_obs[fig3_summary$var == 'grain_yield'], ')'),
    'heading_time' = paste0('Heading time\n(Trials: ', fig3_summary$nenv[fig3_summary$var == 'heading_time'], 
                            ', Obs.: ', fig3_summary$tot_obs[fig3_summary$var == 'heading_time'], ')'),
    'plant_height' = paste0('Plant height\n(Trials: ', fig3_summary$nenv[fig3_summary$var == 'plant_height'], 
                            ', Obs.: ', fig3_summary$tot_obs[fig3_summary$var == 'plant_height'], ')'),
    'test_weight' = paste0('Test weight\n(Trials: ', fig3_summary$nenv[fig3_summary$var == 'test_weight'], 
                           ', Obs.: ', fig3_summary$tot_obs[fig3_summary$var == 'test_weight'], ')')
  ))) +
  scale_fill_gradient(name = 'Number of \ngenotypes',
                      high = '#FF552E', low = '#13294B',
                      na.value = 'white') +
  xlab('Year') + ylab('Location') +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        plot.margin = margin(10, 25, 0, 0),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), 
        panel.grid = element_blank())
ggsave('Figures/Figure3.tif', width = 7, height = 4.5, units = 'in', dpi = 300)

## Figure 4 ----
# Boxplot of average reliabilities within a trial for each trait. The horizontal lines 
# represent the median, and “x” is the mean

blups  |>
  dplyr::select(-c(data)) |>
  group_by(var) |>
  mutate(mrel=format(round(mean(m_rel),2),nsmall=2)) |>
  ggplot(aes(x=var,y=m_rel)) +
  geom_boxplot(fill='#FF552E',color='#13294B',width=0.5) +
  geom_point(aes(y=as.numeric(mrel)),shape=4, size=2, color='#13294B') +
  ylab('Reliability') +
  scale_y_continuous(breaks = seq(0,1,by=0.1),limits=c(0,1),expand=c(0.01,0)) +
  xlab(NULL) +
  scale_x_discrete(labels=c('Grain yield','Heading time',
                            'Plant height','Test weight')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())
ggsave('Figures/Figure4.tif', width = 7, height = 4, units = 'in', dpi=300)

## Figure 5 ----
# Genetic trend of grain yield (A), heading time (B), and test weight (C) in soft red 
# winter wheat from 2001 to 2021, based on data from the University of Illinois wheat 
# breeding program. The orange line represents the genetic trend, while the dotted lines
# delineate the upper and lower boundaries of the 95% confidence interval, which are 
# shaded. Each point corresponds to the estimated average genetic value among the test 
# entries released each year. The genetic trend value, standard error, and significance are 
# also shown

trends <- out_trends |>
  dplyr::select(var, fixeff) |>
  unnest(fixeff) |>
  filter(term%in%c('year_n','gidyr_n','at(check_f, FALSE):gidyr_n')) |>
  mutate(slope=round(solution,2),
         se=round(std_error,2),
         sign=ifelse(pr_chisq<0.001,'***',
                     ifelse(pr_chisq<0.01,'**',
                            ifelse(pr_chisq<0.05,'*',''))))

gt_dat <- out_trends |>
  dplyr::select(var,dat_gen_trend) |>
  unnest(cols = c(dat_gen_trend))

### A ----
slope_gy <- trends |>
  filter(var=='grain_yield'&trend=='gen') |>
  mutate(slope=round(solution,1),
         se=round(std_error,1)) |>
  mutate(label=paste(slope, '±', se, 'kg ha<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

YLD <- gt_dat |> filter(var=='grain_yield') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2001,2021, by = 1)) +
  scale_y_continuous(name = expression('Grain yield (kg ha'^-1*')'),
                     limits = c(4200,5800), expand = c(0,0), breaks = seq(3000, 7000, by = 250)) +
  geom_richtext(data=slope_gy,aes(x=2001,y=5625,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### B ----
slope_tw <- trends |>
  filter(var=='test_weight'&trend=='gen') |>
  mutate(label=paste(slope, '±', se,'g L<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

TW <- gt_dat |> filter(var=='test_weight') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2001, 2021, by = 2)) + # Display only odd years
  scale_y_continuous(name = expression('Test weight (g L'^-1*')'),
                     limits = c(718,775), expand = c(0,0), breaks = seq(720, 770, by = 10)) +
  geom_richtext(data=slope_tw,aes(x=2001,y=770,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### C ----
slope_ht <- trends |>
  filter(var=='heading_time'&trend=='gen') |>
  mutate(label=paste(slope, '±', se,'d yr<sup>-1</sup>', sign)) |>
  glimpse()

HD <- gt_dat |> filter(var=='heading_time') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2001, 2021, by = 2)) + # Display only odd years
  scale_y_continuous(name = expression('Heading time (JD)'),
                     limits = c(129,139), expand = c(0,0), breaks = seq(130, 138, by = 2)) +
  geom_richtext(data=slope_ht,aes(x=2001,y=138,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### Combine A, B, and C ----
combined_plot <- YLD + 
  (HD | TW) + 
  plot_layout(heights = c(2, 1)) + # Adjust the height ratio
  plot_annotation(tag_levels = 'A') # Add A, B, C labels
combined_plot
ggsave('Figures/Figure5.tif', plot = combined_plot, width = 7, height = 7, units = 'in', dpi=300)

# Supplementary Material ----

## Table S1 ----
# Rates of change in genetic, nongenetic, and total phenotypic values for soft red 
# winter wheat from the University of Illinois wheat breeding program from 2005 to 2021 for
# grain yield, plant height, test weight, and heading time

out_trends_17yr |>
  dplyr::select(var, fixeff) |>
  unnest(fixeff) |>
  filter(term%in%c('year_n','gidyr_n','at(check_f, FALSE):gidyr_n')) |>
  mutate(slope=round(solution,2),
         se=round(std_error,2),
         sign=ifelse(pr_chisq<0.001,'***',
                     ifelse(pr_chisq<0.01,'**',
                            ifelse(pr_chisq<0.05,'*','')))) |>
  mutate(value=paste(slope, '±', se, sign)) |>
  dplyr::select(var,trend,value) |>
  pivot_wider(names_from = trend, values_from = value)

## Figure S1 ----
# Genetic trend of grain yield (A), heading time (B), and test weight (C) in soft red 
# winter wheat from 2005 to 2021, based on data from the University of Illinois wheat 
# breeding program. The orange line represents the genetic trend, while the dotted lines
# delineate the upper and lower boundaries of the 95% confidence interval, which are 
# shaded. Each point corresponds to the estimated average genetic value among the test 
# entries released each year. The genetic trend value, standard error, and significance are 
# also shown

trends_17yr <- out_trends_17yr |>
  dplyr::select(var, fixeff) |>
  unnest(fixeff) |>
  filter(term%in%c('year_n','gidyr_n','at(check_f, FALSE):gidyr_n')) |>
  mutate(slope=round(solution,2),
         se=round(std_error,2),
         sign=ifelse(pr_chisq<0.001,'***',
                     ifelse(pr_chisq<0.01,'**',
                            ifelse(pr_chisq<0.05,'*',''))))

gt_dat_17yr <- out_trends_17yr |>
  dplyr::select(var,dat_gen_trend) |>
  unnest(cols = c(dat_gen_trend))

### A ----
slope_gy_17yr <- trends_17yr |>
  filter(var=='grain_yield'&trend=='gen') |>
  mutate(slope=round(solution,1),
         se=round(std_error,1)) |>
  mutate(label=paste(slope, '±', se, 'kg ha<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

YLD_17yr <- gt_dat_17yr |> filter(var=='grain_yield') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2005,2021, by = 1)) +
  scale_y_continuous(name = expression('Grain yield (kg ha'^-1*')'),
                     limits = c(4200,5800), expand = c(0,0), breaks = seq(3000, 7000, by = 250)) +
  geom_richtext(data=slope_gy_17yr,aes(x=2005,y=5625,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### B ----
slope_tw_17yr <- trends_17yr |>
  filter(var=='test_weight'&trend=='gen') |>
  mutate(label=paste(slope, '±', se,'g L<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

TW_17yr <- gt_dat_17yr |> filter(var=='test_weight') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2005, 2021, by = 2)) + # Display only odd years
  scale_y_continuous(name = expression('Test weight (g L'^-1*')'),
                     limits = c(708,775), expand = c(0,0), breaks = seq(710, 770, by = 10)) +
  geom_richtext(data=slope_tw_17yr,aes(x=2005,y=770,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### C ----
slope_ht_17yr <- trends_17yr |>
  filter(var=='heading_time'&trend=='gen') |>
  mutate(label=paste(slope, '±', se,'d yr<sup>-1</sup>', sign)) |>
  glimpse()

HD_17yr <- gt_dat_17yr |> filter(var=='heading_time') |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, linewidth=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= 'Year', breaks = seq(2005, 2021, by = 2)) + # Display only odd years
  scale_y_continuous(name = expression('Heading time (JD)'),
                     limits = c(129,141), expand = c(0,0), breaks = seq(130, 140, by = 2)) +
  geom_richtext(data=slope_ht_17yr,aes(x=2005,y=140,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())

### Combine A, B, and C ----
combined_plot_17yr <- YLD_17yr + 
  (HD_17yr | TW_17yr) + 
  plot_layout(heights = c(2, 1)) + # Adjust the height ratio
  plot_annotation(tag_levels = 'A') # Add A, B, C labels
combined_plot_17yr
ggsave('Figures/FigureS1.tif', plot = combined_plot_17yr, width = 7, height = 7, units = 'in', dpi=300)

# End ----