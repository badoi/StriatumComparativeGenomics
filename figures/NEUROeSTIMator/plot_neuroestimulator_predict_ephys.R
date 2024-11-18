library(tidyverse)
library(ggh4x)
library(here)
library(tidymodels)
library(broom)

DATADIR='data/tidy_data/'

##################################
## 0) prepare the metadata
subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)

##############################################################################
## 1) load the neuroestimulator scores of all cells
save_fn1 = here('data/tidy_data/rdas/human_macaque_mouse_rat_striatum_neuroestimulator.meta.rds')
df_stim = readRDS(save_fn1) %>% filter(Case %in% c('Unaffected', 'Saline')) %>% 
  filter(!Project %in% c('Stanley et al.', 'Chen et al.'), Species == 'Mouse', 
         !integrated_subclusters %in% c('D1-NUDAP'), Region == 'NAcc') %>% 
  mutate(integrated_clusters = factor(integrated_clusters, c('eSPN', 'D1'))) %>% 
  filter(!is.na(integrated_clusters))

threshold = 0.6
df_stim_pb = df_stim %>% 
  group_by(Project, Replicate, Sex, Case, integrated_clusters) %>%
  summarise(num_cells = n(), 
            proportion_active = sum(predicted_activity > threshold) / num_cells, 
            predicted_activity = mean(predicted_activity))

## test the difference in proportion of active cells between D1 and D1/2
df_stim_stats = df_stim_pb %>%
  lm(proportion_active ~ integrated_clusters + Sex + Replicate, data = .) %>% 
  tidy() %>% filter(term == 'integrated_clustersD1') %>% 
  mutate(label = ifelse(p.value < 0.05, '*', 'NS')) %>% mutate(tmp = 'tmp')
df_stim_stats

## add the y position of labels
df_stim_stats= df_stim_pb %>% mutate(tmp = 'tmp') %>% group_by(tmp) %>%
  summarise(y_max = max(proportion_active) * 1.1) %>%
  inner_join(x = df_stim_stats)


############################################################################
### plot the predicted activity of the cells in the neuroestimulator dataset
pp1 = ggplot(df_stim_pb, aes(x = integrated_clusters, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(pch = 21, alpha = 0.8, aes(fill = integrated_clusters, size = num_cells)) +
  geom_text(data = df_stim_stats, vjust = 'inward', size = 2.5, 
            aes(label = label, y = y_max, x = 1.5)) +
  scale_fill_manual(values = c('black', 'gray')) +
  scale_x_discrete(labels = c('D1' = 'D1R', 'eSPN' = 'D1/2R')) +
  theme_minimal(base_size = 5) + scale_size_continuous(range = c(0.5, 2)) +
  labs(x = 'Cluster', y = paste('Proportion Predicted\nActive by snRNA-seq')) +
  theme(legend.position = 'none', axis.title.x = element_blank())

pp1

## save the main plot
pdf(here(FIGDIR, 'pred_ephys_state_neuroestimulator.pdf'), width = 1, height = 1)
pp1
dev.off()



