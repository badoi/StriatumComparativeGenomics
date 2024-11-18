library(tidyverse)
library(ggh4x)
library(here)
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
clusters = c('D1', 'D2', 'eSPN', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,6,4,7)], clusters)

########################################################
## 1) load the neuroestimulator scores of all cells
save_fn1 = here('data/tidy_data/rdas/human_macaque_mouse_rat_striatum_neuroestimulator.meta.rds')
df_stim = readRDS(save_fn1) %>% filter(Case %in% c('Unaffected', 'Saline')) %>% 
  filter(!Project %in% c('Stanley et al.', 'Chen et al.'), 
         integrated_subclusters!= 'D1-NUDAP') %>%
  mutate(Species2 = case_when(Species %in% c('Mouse', 'Rat') ~ 'Mouse & Rat', 
                              Species %in% c('Human', 'Macaque') ~ 'Human & Macaque'), 
         Region2 = case_when(Region %in% c('Caudate', 'Putamen', 'CPu') ~ 
                              'Dorsal Striatum', T ~ 'Nucleus Accumbens'), 
         integrated_clusters = factor(integrated_clusters, clusters))

#############################
## 2) compare patch-matrix
threshold = 0.6
df_stim_pb = df_stim %>% 
  group_by(Project, Species2, Region2, Replicate, Sex, Case,
           integrated_clusters, integrated_subclusters) %>%
  summarise(num_cells = n(), 
            proportion_active = sum(predicted_activity > threshold) / num_cells, 
            predicted_activity = mean(predicted_activity)) %>% 
  filter(integrated_clusters %in% c('D1', 'D2', 'eSPN', 'IC')) %>% 
  mutate(Project = factor(Project, levels = projs) %>% droplevels())

## statistics testing the difference in proportion of active cells between the subclusters
df_stim_stats = df_stim_pb %>% filter(!integrated_clusters %in% c('IC', 'eSPN')) %>% 
  nest(data = -c(Species2, Region2, integrated_clusters)) %>% 
  mutate(mod = map(data, ~ lm( proportion_active ~ integrated_subclusters + 
                                 Replicate + Project , data = .x)),
         tidied = map(mod, tidy)) %>%
  unnest(tidied) %>% select(-mod, -data) %>% 
  filter(grepl('integrated_subclusters', term)) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr')) %>%
  mutate(label = case_when(FDR < 0.001 ~ '***',  FDR < 0.01 ~ '**', 
                           FDR < 0.05 ~ '*', FDR < 0.10 ~ '', T ~ ''))

table(df_stim_stats$label)

## add the y position of labels
df_stim_stats= df_stim_pb %>% group_by(Species2, Region2, integrated_clusters) %>% 
  summarise(y_max = max(proportion_active) * 1.2) %>% inner_join(x = df_stim_stats)

### plot the predicted activity of the cells in the neuroestimulator dataset
pp1 = ggplot(df_stim_pb %>% filter(!integrated_clusters %in% c('IC', 'eSPN')), 
             aes(x = integrated_subclusters, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 21, aes(size =num_cells, fill = integrated_subclusters), alpha = 0.8) +
  geom_text(data = df_stim_stats, aes(label = label, y = y_max, x= 1.5), vjust = 'inward') +
  facet_nested( ~Species2 + Region2 + integrated_clusters, scales = "free", space = 'free_x') +
  scale_fill_manual(values = cols1) +
  theme_classic(base_size = 5) + scale_size_continuous(range = c(0.8, 1.5)) +
  labs(x = 'Cluster', y = paste('Proportion Predicted Active by snRNA-seq')) +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 20, hjust = 1))

## save the main plot
pdf(here(FIGDIR, 'patch_matrix_predicted_activity_neuroestimulator.pdf'), 
    width = 4, height = 1.8)
pp1
dev.off()





#######################################
## 3) compare D1R vs. D2R vs. D1/2R
threshold = 0.6
df_stim_pb = df_stim %>% filter(integrated_clusters %in% c('D1', 'eSPN')) %>%
  group_by(Project, Species2, Region2, Replicate, Sex, Case, integrated_clusters) %>%
  summarise(num_cells = n(), 
            proportion_active = sum(predicted_activity > threshold) / num_cells)

## statistics testing the difference in proportion of active cells between the subclusters
df_stim_stats2 = df_stim_pb %>% filter(integrated_clusters %in% c('D1', 'eSPN')) %>% 
  nest(data = -c(Species2, Region2)) %>% 
  mutate(mod = map(data, ~ lm(proportion_active ~ integrated_clusters + Sex + 
                                Replicate + Project, data = .x)),
         tidied = map(mod, tidy)) %>%
  unnest(tidied) %>% select(-mod, -data) %>% 
  filter(grepl('integrated_clusters', term)) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr')) %>%
  mutate(label = case_when(FDR < 0.001 ~ '***',  FDR < 0.01 ~ '**', 
                           FDR < 0.05 ~ '*', FDR < 0.10 ~ '', T ~ ''), 
         integrated_clusters = ss(term, 'integrated_clusters', 2 ) %>% factor(clusters))

table(df_stim_stats2$label)

## add the y position of labels
df_stim_stats2= df_stim_pb %>% group_by(Species2, Region2, integrated_clusters) %>% 
  summarise(y_max = max(proportion_active) * 1.1) %>% inner_join(x = df_stim_stats2)

### plot the predicted activity of the cells in the neuroestimulator dataset
pp2 = ggplot(df_stim_pb %>% filter(!integrated_clusters %in% c('IC')), 
       aes(x = integrated_clusters, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 21, aes(size =num_cells, fill = integrated_clusters), alpha = 0.8) +
  geom_text(data = df_stim_stats2, vjust = 'inward', x = 1.5,
            aes(label = label, y = y_max)) +
  facet_nested( ~Species2 + Region2, scales = "free", space = 'free') +
  scale_fill_manual(values = cols2) +
  scale_x_discrete(labels = c('D1' = 'D1R', 'D2' = 'D2R', 'eSPN' = 'D1/2R')) +
  theme_classic(base_size = 5) + scale_size_continuous(range = c(0.8, 1.5)) +
  labs(x = 'Cluster', y = paste('Proportion Predicted Active by snRNA-seq')) +
  theme(legend.position = 'none')

## save the main plot
pdf(here(FIGDIR, 'D12R_v_D1R_D2R_predicted_activity_neuroestimulator.pdf'), 
    width = 4, height = 1.8)
pp2
dev.off()

####################################
## 4) export the statistics tables
out_list = list('A' = df_stim_stats, 'B' = df_stim_stats2)

out_fn = here(TABDIR, 'striatum_celltype_predicted_activity_stats.xlsx')
writexl::write_xlsx(out_list, out_fn)