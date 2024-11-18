library(tidyverse)

cols3 = ArchR::paletteDiscrete(c('D1', 'D2', 'eSPN', 'IC'), set = 'paired')
cols4 = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                 c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2-Hybrid', 
                   'D2-Matrix', 'D2-Striosome', 'IC'))

df = readRDS('data/tidy_data/Human_Phan/OUD_Striatum_refined_msn_SeuratObj_N22.neuroestimator.rds') %>% filter(!is.na(neuroestimator)) %>% 
  mutate(activated = ifelse(neuroestimator > 0.5, 'Activated', 'Not activated'), 
         celltypes3 = factor(celltype3, levels =names(cols4)), 
         DSM.IV.OUD = ifelse(DSM.IV.OUD == 'CTL', 'UC', 'OUD'), 
         DSM.IV.OUD = factor(DSM.IV.OUD, levels = c('UC', 'OUD')))

ggplot(df, aes(x = neuroestimator)) + 
  stat_ecdf(geom = 'line', size = 1, aes(color = celltype3, lty =DSM.IV.OUD)) + 
  scale_color_manual(values=cols4) + 
  theme_minimal() + 
  facet_wrap(~Region) + 
  theme(legend.position = 'bottom') + 
  labs(title = 'Neuroestimator integrated subclusters', 
       x = 'Neuroestimator score', 
       y = 'ECDF')


df_sum = df %>% group_by(ID, Case, Region, DSM.IV.OUD, celltype3) %>% 
  summarise(numActivated = sum(activated == 'Activated'), 
            propActivated = numActivated/n()) %>% ungroup()

ggplot(df_sum, aes(x = celltype3, y = propActivated,fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA, alpha = .7) +
  geom_point(pch = 21, size = 3, position = position_jitterdodge()) +
  scale_fill_manual(values=c('black', 'red')) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  labs(x = 'Cell type', y = 'Proportion of MSNs Active')
