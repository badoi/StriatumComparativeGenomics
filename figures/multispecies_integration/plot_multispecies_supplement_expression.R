library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(patchwork)
library(viridis)
library(cluster)
library(future)
library(here)
library(ggh4x)
library(scattermore)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

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

#############################################################
## 1) load the single cell embedding across multiple species
obj_rodents = readRDS(here('data/tidy_data/rdas/integrated_rodent_striatum.supercells.rds'))
obj_primate = readRDS(here('data/tidy_data/rdas/integrated_primate_striatum.supercells.rds'))

obj_rodents$Project = factor(obj_rodents$Project, levels = proj1) %>% droplevels()
obj_primate$Project = factor(obj_primate$Project, levels = proj1) %>% droplevels()

obj_rodents$integrated_clusters = factor(obj_rodents$integrated_clusters, levels = clusters)
obj_primate$integrated_clusters = factor(obj_primate$integrated_clusters, levels = clusters)

obj_rodents$integrated_subclusters = factor(obj_rodents$integrated_subclusters, levels = subclusters)
obj_primate$integrated_subclusters = factor(obj_primate$integrated_subclusters, levels = subclusters)

########################################################################
## 2) plot the expression levels of the genes of interest in primates
# Violin plot - Visualize single cell expression distributions in each cluster
features1 = c('DRD1', 'DRD2', 'DRD3', 'CASZ1', 'STXBP6', 'KCNIP1', 'TSHZ1')
features2 = c('Drd1', 'Drd2', 'Drd3', 'Casz1', 'Stxbp6', 'Kcnip1', 'Tshz1')

df1 = LayerData(obj_primate, layer = 'data')[features1,] %>% t() %>% 
  as.matrix() %>% as.data.frame() %>% bind_cols(obj_primate[[]]) 
df2 = LayerData(obj_rodents, layer = 'data')[features2,] %>% t() %>% 
  as.matrix() %>% as.data.frame() %>% bind_cols(obj_rodents[[]]) %>% 
  rename('DRD1' = 'Drd1', 'DRD2' = 'Drd2', 'DRD3' = 'Drd3', 'CASZ1' = 'Casz1', 
         'STXBP6' = 'Stxbp6', 'KCNIP1' = 'Kcnip1', 'TSHZ1' = 'Tshz1')

df = bind_rows(df1, df2)
df_long = df %>% 
  pivot_longer(cols = all_of(features1), names_to = 'gene', values_to = 'value') %>% 
  mutate(gene = factor(gene, levels = features1))

pdf(here::here(PLOTDIR, 'gene_integrated_clusters.vln.pdf'), 
    width = 8.5, height = 6)
df_long %>% 
  ggplot(aes(x = integrated_subclusters, y = value, fill = integrated_subclusters)) + 
  geom_violin(scale = 'width', width = 0.5, linewidth = 0.25, 
              draw_quantiles = c(0.25, .5, .75)) + 
  scale_fill_manual(values = cols2, guide = 'none') + coord_flip() + 
  scale_x_discrete(limits = rev) + theme_classic(base_size = 5) +
  facet_grid2(gene ~ Project, scales = 'free_x', independent = 'x') 
dev.off()

############################################################
## 3) plot SPN cell type co-expression  - DRD1, DRD2, DRD3
p1 = ggplot(df, aes(x = DRD1, y = DRD2)) +
  geom_point(aes(color = integrated_subclusters, size = size), pch = 20, alpha = .8) +
  facet_wrap(~Project, scales = 'free', nrow = 1) + theme_classic(base_size = 5) + 
  scale_color_manual(values = cols2, guide = 'none') + 
  scale_size_continuous(range = c(0.1, 1), guide = 'none')

p2 = ggplot(df, aes(x = DRD1, y = DRD3)) +
  geom_point(aes(color = integrated_subclusters, size = size), pch = 20, alpha = .8) +
  facet_wrap(~Project, scales = 'free', nrow = 1) + theme_classic(base_size = 5) + 
  scale_color_manual(values = cols2, guide = 'none') + 
  scale_size_continuous(range = c(0.1, 1), guide = 'none')

p3 = ggplot(df, aes(x = DRD2, y = DRD3)) +
  geom_point(aes(color = integrated_subclusters, size = size), pch = 20, alpha = .8) +
  facet_wrap(~Project, scales = 'free', nrow = 1) + theme_classic(base_size = 5) + 
  scale_color_manual(values = cols2, guide = 'none') + 
  scale_size_continuous(range = c(0.1, 1), guide = 'none')

png(here::here(PLOTDIR, 'dopamine_receptor_co-expression_projects.png'), 
    width = 8.5, height = 3, units = 'in', res = 1000)
p1 + p2 + p3 + plot_layout(ncol = 1)
dev.off()

################################################
## 4) plot the per cell QC metrics by project
metadata = bind_rows(obj_rodents[[]], obj_primate[[]]) %>% 
  mutate(nCount_RNA = nCount_RNA / size, nFeature_RNA = nFeature_RNA / size) %>% 
  mutate(Project = factor(Project, levels = proj1))

pdf(here::here(PLOTDIR, 'QC_plots_byProject.vln.pdf'), width = 8.5, height = 2)
p1 = ggplot(metadata, aes(y = nCount_RNA, x = integrated_clusters)) + 
  geom_violin(aes(fill = integrated_clusters), linewidth = 0.25,
              draw_quantiles = c(0.25, .5, .75)) + 
  scale_fill_manual(values = cols1, guide = 'none') +
  facet_grid( ~ Project) + theme_classic(base_size = 5) + 
  ylab('# Detected UMI per cell')+ xlab('') + scale_y_log10()

p2 = ggplot(metadata, aes(y = nFeature_RNA, x = integrated_clusters)) + 
  geom_violin(aes(fill = integrated_clusters),  linewidth = 0.25,
              draw_quantiles = c(0.25, .5, .75)) + 
  scale_fill_manual(values = cols1, guide = 'none') +
  facet_grid( ~ Project) + theme_classic(base_size = 5) +
  ylab('# Detected Genes per cell') + xlab('') + scale_y_log10() 

p1 + p2 + plot_layout(ncol = 1) & theme(plot.tag = element_text(size = 5))
dev.off()

