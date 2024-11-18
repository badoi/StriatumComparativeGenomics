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
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
TABLDIR = 'figures/mouse_tshz1_pdyn/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

##################################
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

#############################################################
## 1) load the single cell embedding across multiple species
supercells_fn = 
  here('data/tidy_data/Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.supercells.rds')
obj_mouse = readRDS(supercells_fn) %>% subset(Region %in% 'NAcc') 

##################################################################
## 2) plot the integrated subclusters labels for the mouse NAcc
mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', plot.title = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_mouse, feature1 = "Pdyn", feature2 = "Tshz1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'mouse_oprm1_pdyn_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_mouse, feature1 = "Pdyn", feature2 = "Oprm1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'mouse_oprm1_tshz1_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_mouse, feature1 = "Tshz1", feature2 = "Oprm1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()


###############################################################################
## 3) count division of integrated subclusters by Tshz1+ or Pdyn+ 
thresh1 = 1.5; thresh2 = 1.0
keep_groups = c('Tshz1+, Pdyn-', 'Tshz1-, Pdyn+', 'Tshz1-, Pdyn-')

df = LayerData(obj_mouse, 'data')[c('Tshz1', 'Pdyn', 'Oprm1'),] %>%as.matrix() %>% 
  t() %>% as.data.frame() %>% bind_cols(x = obj_mouse[[]]) %>% 
  mutate(Tshz1_group = ifelse(Tshz1 > thresh1, 'Tshz1+', 'Tshz1-'),
         Tshz1_group = factor(Tshz1_group, levels = c('Tshz1+', 'Tshz1-')),
         Pdyn_group = ifelse(Pdyn > thresh1, 'Pdyn+', 'Pdyn-'),
         Pdyn_group = factor(Pdyn_group, levels = c('Pdyn-', 'Pdyn+')),
         Oprm1_group = ifelse(Oprm1 > thresh2, 'Oprm1+', 'Oprm1-'),
         Oprm1_group = factor(Oprm1_group, levels = c('Oprm1-', 'Oprm1+')), 
         group = paste(Tshz1_group, Pdyn_group, sep = ', ') %>% factor())

df_group = df %>% filter(group %in% keep_groups) %>% 
  group_by(group, integrated_subclusters) %>% 
  summarise(n = n()) %>%
  mutate(integrated_subclusters = factor(integrated_subclusters, subclusters)) %>% 
  arrange(desc(integrated_subclusters)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_integrated_subclusters.pie.pdf'), 
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


df_group2 = df %>% group_by(group) %>% summarise(n = n()) %>%
  arrange(desc(group)) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.01, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_overlap.pie.pdf'), height = 1.5, width = 1.2)
ggplot(df_group2, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = group)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos), size = 2) +
  scale_fill_manual(values = c('gray', '#1f78b4', '#b856d7',  'black'), name ='') + 
  theme_classic(base_size = 5) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position = 'bottom', 
        legend.key.width = unit(.5, "lines"), 
        legend.key.height = unit(.5, "lines"))
dev.off()

##############################################################################
## 4) count division of integrated subclusters by Tshz1+ or Pdyn+ and Oprm1
df_group2 = df %>% group_by(group, Oprm1_group) %>% 
  summarise(n = n()) %>%
  arrange(desc(Oprm1_group)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_oprm1_proportion.pie.pdf'), 
    height = 3, width = 1.2)
ggplot(df_group2, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = Oprm1_group)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos), size = 2) +
  scale_fill_manual(values = c('gray', 'darkgreen'), name = '') + 
  theme_classic(base_size = 5) + facet_wrap( ~group, ncol =1 ) + 
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position = 'bottom')
dev.off()


##############################################################################
## 5) quantify the per-cell expression of Oprm1 across Tshz1+ and Pdyn+ cells
levels(df$group) = str_wrap(levels(df$group), 5)

pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_oprm1_expression.vln.pdf'), 
    height = 1, width = 1.2)
ggplot(df, aes( x=group, y= Oprm1)) + 
  geom_violin(aes(fill = group), draw_quantiles = c(.25, .5, .75)) +
  scale_fill_manual(values = c('#b856d7', '#1f78b4', 'grey'), name = '') + 
  theme_classic(base_size = 5) + 
  theme(legend.position = 'none')
dev.off()

df_test = lmer(Oprm1 ~ group + Sex + (1|Replicate), data = df) %>% tidy() %>% 
  filter(effect == 'fixed', str_detect(term, 'group')) %>% 
  mutate(term = str_replace(term, 'group', 'Tshz1+, Pdyn- vs. ') %>% str_replace_all('\n', ' '),, 
         estimate = -estimate, statistic = -statistic, 
         p.bonf = p.adjust(p.value, 'bonferroni')) %>% 
  as.data.frame() %>% 
  select(-effect, -group)

#                              term  estimate  std.error statistic   df
# 1 Tshz1+, Pdyn- vs. Tshz1-, Pdyn+ 0.3575483 0.02312898  15.45888 5239
# 2 Tshz1+, Pdyn- vs. Tshz1-, Pdyn- 0.6519039 0.01922633  33.90683 5239
#         p.value        p.bonf
# 1  9.465988e-53  1.893198e-52
# 2 5.155129e-228 1.031026e-227
