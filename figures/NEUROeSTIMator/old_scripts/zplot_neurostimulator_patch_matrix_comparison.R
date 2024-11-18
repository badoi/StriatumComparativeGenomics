library(tidyverse)
library(ggh4x)
library(here)
library(outliers)
library(patchwork)
library(tidymodels)
library(broom)

DATADIR='data/tidy_data/'
FIGDIR='figures/NEUROeSTIMator/plots/'
TABDIR='figures/NEUROeSTIMator/tables/'
dir.create(here(FIGDIR), showWarnings = FALSE, recursive = TRUE)
dir.create(here(TABDIR), showWarnings = FALSE, recursive = TRUE)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

##############################################################################
## 0) prepare the metadata
projs = c('Gayden Kozel et al.','Siletti et al.', 'Tran et al.', 'Phan et al.', 
  'Li et al.', 'Corces et al.', 'Chiou et al.', 'He, Kleyman et al.', 
  'Phillips et al.', 'Savell et al.', 'Zeng et al.', 'Chen et al.',
  'Saunders et al.',  'Stanley et al.',  'Zu et al.')

subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)

##############################################################################
## 1) load the neuroestimulator scores of all cells
save_fn1 = here('data/tidy_data/rdas/human_macaque_mouse_rat_striatum_neuroestimulator.meta.rds')
df_stim = readRDS(save_fn1) %>% filter(Case %in% c('Unaffected', 'Saline')) %>% 
  filter(!Project %in% c('Stanley et al.')) %>%
  mutate(Project = as.character(Project), 
         Project = ifelse(Project == 'Savell et al.', 'Phillips et al.', Project),
         Project = factor(Project, levels = projs) %>% droplevels(),
         Species2 = case_when(Species %in% c('Mouse', 'Rat') ~ 'Mouse & Rat', 
                              Species %in% c('Human', 'Macaque') ~ 'Human & Macaque'), 
         Region2 = case_when(Region %in% c('Caudate', 'Putamen', 'CPu') ~ 
                              'Dorsal Striatum', T ~ 'Nucleus Accumbens'))

