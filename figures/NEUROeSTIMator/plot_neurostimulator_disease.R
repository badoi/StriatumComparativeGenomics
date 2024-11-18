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
df_stim = readRDS(save_fn1) 

## calculate the proportion of active cells from single cell predicted activity
threshold = 0.6
df_stim_pb = df_stim %>% 
  group_by(Project, Species, Region, Replicate, Sex, Case,
           integrated_clusters, integrated_subclusters) %>%
  summarise(num_cells = n(), 
            proportion_active = sum(predicted_activity > threshold) / num_cells, 
            predicted_activity = mean(predicted_activity)) %>% 
  filter(integrated_clusters %in% c('D1', 'D2', 'eSPN', 'IC')) %>% 
  mutate(Project = factor(Project, levels = projs) %>% droplevels(), 
         integrated_clusters = factor(integrated_clusters, levels = clusters))

## statistics testing the difference in proportion of active cells between the subclusters
df_stim_cocaine = df_stim_pb %>%
  filter(Project %in% c('Phillips et al.', 'Savell et al.')) %>% 
  mutate(Case = factor(Case, c('Saline', 'Cocaine-Acute', 'Cocaine-Repeated'))) %>% 
  group_by(integrated_clusters) %>% nest() %>% 
  mutate(mod = map(data, ~
                     if(length(unique(.x$integrated_subclusters)) > 1){
                       lm(proportion_active ~ Case + Sex+ integrated_subclusters, data = .x)
                     } else {
                       lm(proportion_active ~ Case + Sex, data = .x)
                     }),
         tidied = map(mod, tidy)) %>%
  unnest(tidied) %>% select(-mod, -data) %>%
  filter(grepl('Case', term)) %>% 
  mutate(label = case_when(p.value < 0.001 ~ '***',  p.value < 0.01 ~ '**', 
                           p.value < 0.05 ~ '*', p.value < 0.10 ~ '.', 
                           T ~ 'ns'), 
         Case = str_remove(term, 'Case'))

table(df_stim_cocaine$label)

## add the y position of labels
df_stim_cocaine = df_stim_pb %>% 
  filter(Project  %in% c('Phillips et al.', 'Savell et al.')) %>% 
  group_by(integrated_clusters) %>% 
  summarise(y_max = max(proportion_active) * 1.2) %>%
  inner_join(x = df_stim_cocaine)

## make the plot of activity by cocaine exposure
pp3 = df_stim_pb %>% filter(Project  %in% c('Phillips et al.', 'Savell et al.')) %>% 
  mutate(Case = factor(Case, levels = c('Saline', 'Cocaine-Acute', 'Cocaine-Repeated'))) %>%
  ggplot(aes(x = Case, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(dodge.width=.5), pch = 24, 
             aes(fill = integrated_subclusters, size = num_cells), alpha = 0.8) +
  geom_text(data = df_stim_cocaine, vjust = 'inward', size = 2,
            aes(label = label, y = y_max, x = Case)) +
  facet_nested_wrap( ~ "Savell et al. & Phillips et al., Rat, Nucleus Accumbens" + 
                       integrated_clusters, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = cols1) + scale_size_continuous(range = c(1, 1.8)) +
  theme_minimal(base_size = 5) + scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Condition', y = paste('Proportion Predicted\nActive by snRNA-seq')) +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 30, hjust = 1))

##########################################################################################
## statistics testing the difference in proportion of active cells between the subclusters
df_stim_oud = df_stim_pb %>% filter(Project == 'Phan et al.') %>%
  filter(!integrated_subclusters %in% c('D1-NUDAP')) %>%
  lm(proportion_active ~ integrated_clusters + Sex + Region + num_cells + 
       integrated_clusters:Case, data = .) %>% 
  tidy() %>% filter(grepl(':Case', term)) %>% 
  mutate(label = case_when(p.value < 0.001 ~ '***',  p.value < 0.01 ~ '**', 
                           p.value < 0.05 ~ '*', p.value < 0.10 ~ '.',
                           T ~ 'ns'), 
         Case = term %>% ss(':', 2) %>% str_remove('Case'), 
         integrated_clusters = term %>% ss(':') %>% str_remove('integrated_clusters'), 
         integrated_clusters = factor(integrated_clusters, levels = clusters))

## add the y position of labels
df_stim_oud = df_stim_pb %>% filter(Project == 'Phan et al.') %>% 
  group_by(integrated_clusters) %>% 
  summarise(y_max = max(proportion_active) * 0.8) %>% inner_join(x = df_stim_oud)

table(df_stim_oud$label)

pp4 = ggplot(df_stim_pb2 %>% filter(!integrated_subclusters %in% c('D1-NUDAP')),
             aes(x = Case, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0, NA)) +
  geom_point(position=position_jitterdodge(dodge.width=.5), alpha = 0.8, pch = 24,
             aes(fill = integrated_subclusters, size = num_cells, 
                 shape = Region)) +
  geom_text(data = df_stim_oud, vjust = 'inward', size = 2,
            aes(label = label, y = y_max, x = 1.5)) +
  facet_nested_wrap( ~ 'Phan et al., Human, Dorsal Striatum' + integrated_clusters,
                     scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cols1) + 
  theme_minimal(base_size = 5) + scale_size_continuous(range = c(1, 1.8)) +
  labs(x = 'Condition', y = paste('Proportion Predicted\nActive by snRNA-seq')) +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 30, hjust = 1))


## save the main plot
pdf(here(FIGDIR, 'disease_state_neuroestimulator.pdf'), width = 4, height = 2)
pp3 + pp4 + plot_layout(nrow = 1, widths = c(1.5, 1))
dev.off()



#####################
# export the statistics tables
out_list = list('C' = df_stim_cocaine, 'D' = df_stim_oud)

out_fn = here(TABDIR, 'disease_state_predicted_activity_stats.xlsx')
writexl::write_xlsx(out_list, out_fn)

