library(tidyverse)

cols3 = ArchR::paletteDiscrete(c('D1', 'D2', 'eSPN', 'IC'), set = 'paired')
cols4 = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                 c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2-Hybrid', 
                   'D2-Matrix', 'D2-Striosome', 'IC'))

df = readRDS('data/tidy_data/Rat_Phillips/Phillips2023_snRNA_filtered_SeuratObj_N8.neuroestimator.rds') %>%
  filter(!is.na(neuroestimator), grepl('MSN', Combo_CellType)) %>% 
  mutate(activated = ifelse(neuroestimator > 0.5, 'Activated', 'Not activated'), 
         Stim = ifelse(Stim == 'Coc' & Dataset =='Acute', 'Cocaine-Acute', 
                ifelse(Stim == 'Coc' & Dataset =='Repeated', 
                       'Cocaine-Repeated', 'Saline')),
         Stim = factor(Stim, c('Saline', 'Cocaine-Acute', 'Cocaine-Repeated')), 
         Combo_CellType = str_replace_all(Combo_CellType, '-[12]', ''))

ggplot(df, aes(x = neuroestimator)) + 
  stat_ecdf(geom = 'line', size = 1, aes(color = Combo_CellType, lty =Stim)) + 
  scale_color_brewer(palette = 'Paired') + 
  theme_minimal() + 
  facet_wrap(~Dataset) + 
  theme(legend.position = 'bottom') + 
  labs(title = 'Neuroestimator integrated subclusters', 
       x = 'Neuroestimator score', 
       y = 'ECDF')


df_sum = df %>% group_by(orig.ident, Sex, Stim, Dataset, Combo_CellType) %>% 
  summarise(numActivated = sum(activated == 'Activated'), 
            propActivated = numActivated/n()) %>% ungroup()

ggplot(df_sum, aes(x = Combo_CellType, y = propActivated,fill = Stim)) + 
  geom_boxplot(outlier.shape = NA, alpha = .7) +
  geom_point(pch = 21, size = 3, position = position_jitterdodge()) +
  scale_fill_manual(values=c('black', 'red', 'blue')) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,NA)) +
  theme(legend.position = 'bottom') + 
  labs(x = 'Cell type', y = 'Proportion of MSNs Active')
