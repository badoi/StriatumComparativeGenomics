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
df_stim = readRDS(save_fn1) %>% 
  mutate(Project = as.character(Project), 
         Project = ifelse(Project == 'Savell et al.', 'Phillips et al.', Project),
         Project = factor(Project, levels = projs) %>% droplevels(),
         Species2 = case_when(Species %in% c('Mouse', 'Rat') ~ 'Rodent', 
                              Species %in% c('Human', 'Macaque') ~ 'Primate'), 
         Region2 = case_when(Region %in% c('Caudate', 'Putamen', 'CPu') ~ 
                              'Dorsal Striatum', T ~ 'Nucleus Accumbens'))

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
df_stim_stats = df_stim_pb %>% 
  filter(integrated_clusters != 'IC', Case %in% c('Unaffected', 'Saline')) %>% 
  group_by(Project, Species2, Region2, integrated_clusters) %>% 
  filter(length(unique(integrated_subclusters))== 2) %>% ungroup() %>% 
  nest(data = -c(Project, Species2, Region2, integrated_clusters)) %>% 
  mutate(mod = map(data, ~ if(length(unique(.x$Replicate)) > 1){
    lm(proportion_active ~ integrated_subclusters + Replicate, data = .x)
    } else {
      lm(proportion_active ~ integrated_subclusters, data = .x)
    }),
         tidied = map(mod, tidy)) %>%
  unnest(tidied) %>% select(-mod, -data) %>% 
  filter(grepl('integrated_subclusters', term)) %>% 
  group_by(Project) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr')) %>% ungroup() %>% 
  mutate(label = case_when(FDR < 0.001 ~ '***',  FDR < 0.01 ~ '**', 
                           FDR < 0.05 ~ '*', FDR < 0.10 ~ '.', T ~ ''), 
         Project = factor(Project, levels = projs) %>% droplevels())

table(df_stim_stats$label)

## add the y position of labels
df_stim_stats= df_stim_pb %>% 
  group_by(Project, Species2, Region2, integrated_clusters) %>% 
  summarise(y_max = max(proportion_active) * 1.25) %>%
  inner_join(x = df_stim_stats)


############################################################################
### plot the predicted activity of the cells in the neuroestimulator dataset
pp1 = ggplot(df_stim_pb %>% filter(Species2 == 'Rodent', 
                                   Case %in% c('Unaffected', 'Saline')), 
       aes(x = integrated_clusters, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(dodge.width=.9), pch = 21, 
             aes(size =num_cells, fill = integrated_subclusters), alpha = 0.8) +
  geom_text(data = df_stim_stats %>% filter(Species2 == 'Rodent'),
            aes(label = label, y = y_max), vjust = 'inward') +
  facet_nested_wrap( ~ Region2 + Project, scales = "free", nrow = 1) +
  scale_fill_manual(values = cols1) +
  theme_bw(base_size = 5) + scale_size_continuous(range = c(1, 3)) +
  labs(x = 'Cluster', y = paste('Proportion Active')) +
  theme(legend.position = 'none')

pp2 = ggplot(df_stim_pb %>% filter(Species2 == 'Primate', Case == 'Unaffected'), 
       aes(x = integrated_clusters, y = proportion_active)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(dodge.width=.9), pch = 21, 
             aes(size =num_cells, fill = integrated_subclusters), alpha = 0.8) +
  geom_text(data = df_stim_stats %>% filter(Species2 == 'Primate'), vjust = 'inward',
            aes(label = label, y = y_max, x = integrated_clusters)) +
  facet_nested_wrap( ~ Region2 + Project, scales = "free", nrow = 1) +
  scale_fill_manual(values = cols1) +
  theme_bw(base_size = 5) + scale_size_continuous(range = c(1, 3)) +
  labs(x = 'Cluster', y = paste('Proportion Active')) +
  theme(legend.position = 'bottom',
        legend.margin=margin(rep(0, 4)),
        legend.box.margin=margin(rep(-5, 4))) +
  guides(fill = guide_legend(nrow = 1)) 


## save the main plot
pdf(here(FIGDIR, 'unaffected_state_neuroestimulator.pdf'), width = 7, height = 3)
pp1+ pp2 + plot_layout(ncol = 1, heights = c(1, 1.1))
dev.off()

summary(ggplot_build(pp2)$data[[2]]$size)
summary(df_stim_pb$num_cells)


#####################
# export the statistics tables
out_list = list('A' = df_stim_stats %>% filter(Species2 %in% c('Rodent')), 
                'B' = df_stim_stats %>% filter(Species2 %in% c('Primate')))

out_fn = here(TABDIR, 'striatum_celltype_predicted_activity_stats.xlsx')
writexl::write_xlsx(out_list, out_fn)
