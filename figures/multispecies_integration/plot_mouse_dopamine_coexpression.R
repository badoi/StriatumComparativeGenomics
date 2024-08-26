library(tidyverse)
library(Seurat)
library(data.table)
library(patchwork)
library(cluster)
library(future)
library(here)
library(AUCell)
library(lme4)
library(lmerTest)
library(broom.mixed)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
TABLDIR = 'figures/multispecies_integration/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

############################
## 0) prepare the metadata
proj1 = c('Gayden Kozel et al.','Siletti et al.', 'Tran et al.', 'Phan et al.', 
          'Li et al.', 'Corces et al.', 'Chiou et al.', 'He, Kleyman et al.', 
          'Phillips et al.', 'Savell et al.', 'Zeng et al.', 'Chen et al.',
          'Saunders et al.',  'Stanley et al.',  'Zu et al.')

clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)
subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

###################################################
## 1) load the single cell data for the mouse NAcc
supercells_fn = 
  here('data/tidy_data/rdas/integrated_rodent_striatum.supercells.rds')
obj_super = readRDS(supercells_fn) 

##################################################################
## 2) plot the integrated subclusters labels for the mouse NAcc
mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', plot.title = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

pdf(here(PLOTDIR, 'mouse_drd1_drd2_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "Drd2", feature2 = "Drd1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'mouse_drd3_drd2_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "Drd2", feature2 = "Drd3", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'mouse_drd3_drd1_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "Drd1", feature2 = "Drd3", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()


###############################################################################
## 3) count division of integrated subclusters by drd1+ or drd2+ 
thresh1 = 1.5; thresh2 = 1.0
keep_groups = c('Drd1+, drd2-', 'Drd1-, drd2+', 'Drd1-, drd2-')

df = LayerData(obj_super, 'data')[c('Drd1', 'Drd2', 'Drd3'),] %>%as.matrix() %>% 
  t() %>% as.data.frame() %>% bind_cols(x = obj_super[[]]) %>% 
  mutate(drd1_group = ifelse(drd1 > thresh1, 'Drd1+', 'Drd1-'),
         drd1_group = factor(drd1_group, levels = c('Drd1+', 'Drd1-')),
         drd2_group = ifelse(drd2 > thresh1, 'Drd2+', 'Drd2-'),
         drd2_group = factor(drd2_group, levels = c('Drd2-', 'Drd2+')),
         drd3_group = ifelse(drd3 > thresh2, 'Drd3+', 'Drd3-'),
         drd3_group = factor(drd3_group, levels = c('Drd3-', 'Drd3+')), 
         group = paste(drd1_group, drd2_group, sep = ', '), 
         group = factor(group, keep_groups)) %>% 
  filter(group %in% keep_groups)

df_group = df %>% group_by(group, integrated_subclusters) %>% 
  summarise(n = n()) %>%
  mutate(integrated_subclusters = factor(integrated_subclusters, subclusters)) %>% 
  arrange(desc(integrated_subclusters)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'mouse_drd1_drd2_integrated_subclusters.pie.pdf'), 
    height = 3, width = 1)
ggplot(df_group, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = integrated_subclusters)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos), size = 2) +
  scale_fill_manual(values = cols2, guide = 'none') + 
  theme_classic(base_size = 5) + facet_wrap( ~group, ncol =1 ) + 
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank())
dev.off()


##############################################################################
## 4) count division of integrated subclusters by drd1+ or drd2+ and drd3
df_group2 = df %>% group_by(group, drd3_group) %>% 
  summarise(n = n()) %>%
  arrange(desc(drd3_group)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'mouse_drd1_drd2_drd3_proportion.pie.pdf'), 
    height = 3, width = 1.2)
ggplot(df_group2, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = drd3_group)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos), size = 2) +
  scale_fill_manual(values = c('gray', 'darkgreen'), name = '') + 
  theme_classic(base_size = 5) + facet_wrap( ~group, ncol =1 ) + 
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position = 'bottom')
dev.off()


##############################################################################
## 5) quantify the per-cell expression of drd3 across drd1+ and drd2+ cells
levels(df$group) = str_wrap(levels(df$group), 5)

pdf(here(PLOTDIR, 'mouse_drd1_drd2_drd3_expression.vln.pdf'), 
    height = 1, width = 1.2)
ggplot(df, aes( x=group, y= drd3)) + 
  geom_violin(aes(fill = group), draw_quantiles = c(.25, .5, .75)) +
  scale_fill_manual(values = c('#33a02c', '#1f78b4', 'grey'), name = '') + 
  theme_classic(base_size = 5) + 
  theme(legend.position = 'none')
dev.off()

df_test = lmer(drd3 ~ group + Sex + (1|Replicate), data = df) %>% tidy() %>% 
  filter(effect == 'fixed', str_detect(term, 'group')) %>% 
  mutate(term = str_replace(term, 'group', 'Drd1+, drd2- vs. ') %>% str_replace_all('\n', ' '),, 
         estimate = -estimate, statistic = -statistic, 
         p.bonf = p.adjust(p.value, 'bonferroni')) %>% 
  as.data.frame() %>% 
  select(-effect, -group)

#                              term  estimate  std.error statistic   df
# 1 drd1+, drd2- vs. drd1-, drd2+ 0.3575483 0.02312898  15.45888 5239
# 2 drd1+, drd2- vs. drd1-, drd2- 0.6519039 0.01922633  33.90683 5239
#         p.value        p.bonf
# 1  9.465988e-53  1.893198e-52
# 2 5.155129e-228 1.031026e-227
