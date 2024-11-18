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
obj_super = here('data/tidy_data/rdas/integrated_primate_striatum.supercells.rds') %>% 
  readRDS() %>% subset(Region %in% 'NAcc') %>% subset(Species == 'Human') 


##################################################################
## 2) plot the integrated subclusters labels for the human NAcc
mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', plot.title = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

pdf(here(PLOTDIR, 'human_tshz1_pdyn_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "PDYN", feature2 = "TSHZ1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'human_oprm1_pdyn_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "PDYN", feature2 = "OPRM1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()

pdf(here(PLOTDIR, 'human_oprm1_tshz1_co-expression.scatter.pdf'), 
    height = 1.5, width = 1.5)
FeatureScatter(obj_super, feature1 = "TSHZ1", feature2 = "OPRM1", 
               group.by = 'integrated_subclusters', cols = cols2, 
               raster = F, pt.size = .25) & mytheme
dev.off()


###############################################################################
## 3) count division of integrated subclusters by TSHZ1+ or PDYN+ 
## use the co-expression scatter plot to determine the thresholds
thresh1 = 1; thresh2 = .5  
keep_groups = c('TSHZ1+, PDYN-', 'TSHZ1-, PDYN+', 'TSHZ1-, PDYN-')

df = LayerData(obj_super, 'data')[c('TSHZ1', 'PDYN', 'OPRM1'),] %>%as.matrix() %>% 
  t() %>% as.data.frame() %>% bind_cols(x = obj_super[[]]) %>% 
  mutate(TSHZ1_group = ifelse(TSHZ1 > thresh1, 'TSHZ1+', 'TSHZ1-'),
         TSHZ1_group = factor(TSHZ1_group, levels = c('TSHZ1+', 'TSHZ1-')),
         PDYN_group = ifelse(PDYN > thresh2, 'PDYN+', 'PDYN-'),
         PDYN_group = factor(PDYN_group, levels = c('PDYN-', 'PDYN+')),
         OPRM1_group = ifelse(OPRM1 > thresh1, 'OPRM1+', 'OPRM1-'),
         OPRM1_group = factor(OPRM1_group, levels = c('OPRM1-', 'OPRM1+')), 
         group = paste(TSHZ1_group, PDYN_group, sep = ', ') %>% factor())

df_group = df %>% filter(group %in% keep_groups) %>% 
  group_by(group, integrated_subclusters) %>% 
  summarise(n = n()) %>%
  mutate(integrated_subclusters = factor(integrated_subclusters, subclusters)) %>% 
  arrange(desc(integrated_subclusters)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'human_tshz1_pdyn_integrated_subclusters.pie.pdf'), 
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

pdf(here(PLOTDIR, 'human_tshz1_pdyn_overlap.pie.pdf'), height = 1.5, width = 1.2)
ggplot(df_group2, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = group)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos, x = 1.65), size = 2) +
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
## 4) count division of integrated subclusters by TSHZ1+ or PDYN+ and OPRM1
df_group2 = df %>% group_by(group, OPRM1_group) %>% 
  summarise(n = n()) %>%
  arrange(desc(OPRM1_group)) %>% group_by(group) %>% 
  mutate(pc = n/sum(n), ypos = cumsum(pc)- 0.5 * pc) %>% ungroup() %>% 
  mutate(label = paste0(round(pc*100, 0), '%'),
         label = ifelse(pc < 0.05, '', label)) %>% ungroup()

pdf(here(PLOTDIR, 'human_tshz1_pdyn_oprm1_proportion.pie.pdf'), 
    height = 3, width = 1.2)
ggplot(df_group2, aes( x='', y= pc)) + 
  geom_bar(stat = 'identity', aes(fill = OPRM1_group)) +
  coord_polar("y", start=0) + 
  geom_text(aes(label = label, y = ypos), size = 2) +
  scale_fill_manual(values = c('gray', 'darkgreen'), name = '') + 
  theme_classic(base_size = 5) + facet_wrap( ~group, ncol =1 ) + 
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), 
        axis.ticks=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position = 'bottom')
dev.off()


##############################################################################
## 5) quantify the per-cell expression of OPRM1 across TSHZ1+ and PDYN+ cells
levels(df$group) = str_wrap(levels(df$group), 5)

pdf(here(PLOTDIR, 'human_tshz1_pdyn_oprm1_expression.vln.pdf'), 
    height = 1, width = 1.2)
ggplot(df, aes( x=group, y= OPRM1)) + 
  geom_violin(aes(fill = group), draw_quantiles = c(.25, .5, .75)) +
  scale_fill_manual(values = c('#b856d7', '#1f78b4', 'grey'), name = '') + 
  theme_classic(base_size = 5) + 
  theme(legend.position = 'none')
dev.off()

df_test = lmer(OPRM1 ~ group + Sex + Project  + (1|Replicate), 
               data = df) %>% tidy() %>% 
  filter(effect == 'fixed', str_detect(term, 'group')) %>% 
  mutate(term = str_replace(term, 'group', 'TSHZ1+, PDYN- vs. ') %>% str_replace_all('\n', ' '),, 
         estimate = -estimate, statistic = -statistic, 
         p.bonf = p.adjust(p.value, 'bonferroni')) %>% 
  as.data.frame() %>% 
  select(-effect, -group)

#                              term  estimate  std.error statistic       df
# 1 TSHZ1+, PDYN- vs. TSHZ1-, PDYN+ 0.6261133 0.03295490  18.99910 3068.574
# 2 TSHZ1+, PDYN- vs. TSHZ1-, PDYN- 0.6067532 0.02419672  25.07584 3069.181
#         p.value        p.bonf
# 1  3.477181e-76  6.954362e-76
# 2 2.136624e-126 4.273249e-126


