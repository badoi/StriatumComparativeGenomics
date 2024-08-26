library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(viridis)
library(cluster)
library(ggh4x)
library(here)
library(tidymodels)
library(broom)
library(dplyr)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
TABLDIR = 'figures/multispecies_integration/tables'

dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

#################################################################
## 1) calculate the mixing scores from the multi-species datasets
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.RNA.h5Seurat')
obj = LoadH5Seurat(save_fn)
DefaultAssay(obj) = 'integrated'

#################################################
## 2) calculate the proportion of each cell type
table(obj$Project, obj$Region)

df_byDonor = obj[[]] %>% filter(integrated_clusters != 'Other') %>% 
  mutate(Region = ifelse(Region == 'OT', 'NAcc', Region), 
         Region = ifelse(Region == 'NAcc', 'NAcc/OT', Region)) %>% 
  group_by(Project, Species, Region, Replicate, integrated_clusters) %>% 
  dplyr::summarise(num = n()) %>% 
  group_by(Project, Species, Region, Replicate) %>% 
  dplyr::mutate(total = sum(num)) %>% 
  group_by(Project, Species, Region, Replicate, integrated_clusters) %>% 
  dplyr::mutate(prop = num/total * 100) %>% ungroup()
  
df_group = df_byDonor %>% 
  group_by(Species, Region, integrated_clusters) %>% 
  dplyr::summarise(prop_mean = mean(prop), 
                   prop_sem = sd(prop)/sqrt(n())) %>% ungroup() %>% 
  mutate(Region = factor(Region, c( 'Caudate', 'Putamen', 'CPu','NAcc/OT')), 
         Species = factor(Species, c("Human", 'Macaque', 'Mouse', 'Rat')))

stats_df = df_byDonor %>% nest(data = -c(integrated_clusters)) %>% 
  mutate(
    mod = map(data, ~ aov(prop ~ Region + Species, data = .x)), 
    mod2 = map(mod, TukeyHSD), 
    tidy = map(mod2, tidy)
  ) %>% dplyr::select(-c(mod, mod2, data)) %>% 
  unnest(tidy)

out_fn = here('figures/multispecies_integration/tables', 'cell_type_proportion_by_species.xlsx')
writexl::write_xlsx(list(prop = df_group, stats = stats_df), out_fn)


######################################
## 3) plot the cell type proportions
pdf(here::here(PLOTDIR,'cross_species_celltype_proportions.bar.pdf'), 
    width = 7, height = 1.5)

ggplot(df_group, aes(x = Region, y = prop_mean)) + 
  geom_bar(stat = 'identity', aes(fill = Region)) + 
  geom_errorbar(aes(ymin=prop_mean-prop_sem, ymax=prop_mean+prop_sem), width=.4) +
  scale_fill_viridis(discrete=TRUE) + 
  facet_nested( ~ integrated_clusters + Species, scales = 'free_x', space = 'free_x') + 
  theme_bw(base_size = 6) + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  labs(y = 'Percent of Neurons')

dev.off()








