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

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

#############################################################
## 1) load the single cell embedding across multiple species
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.ARC.h5Seurat')
obj = LoadH5Seurat(save_fn)
Idents(obj) = 'Project'

DefaultAssay(obj) = 'RNA'
obj = obj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  
table(obj$celltype.orig, obj$Project)

###############################################################
## 2) make visualization plots of the multi-species integration
DefaultAssay(obj) = 'integrated'
cols = ArchR::paletteDiscrete(obj$celltype.orig)

p1 = DimPlot(obj, reduction = "umap", split.by = "Project", cols = cols, 
        group.by = 'celltype.orig', label = T, label.size = 1.5) +
  guides(colour = guide_legend(nrow = 4, override.aes = list(size = 1)), byrow = TRUE) +
  theme(legend.position = 'bottom', legend.text=element_text(size=6)) + 
  ggtitle("Multi-species striatum datasets: original publication labels") &
  FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
  theme(strip.text.x = element_text(size = 5), 
        legend.margin=margin(0,-20,0,0), legend.box.spacing = unit(0, "pt"),
        legend.spacing.x = unit(-1.5, 'mm'), legend.spacing.y = unit(-2, 'cm')
  )

## plot the main clusters
cols3 = ArchR::paletteDiscrete(obj$integrated_clusters, set = 'paired')
p2 = DimPlot(obj, reduction = "umap", split.by = "Species", cols = cols3,
        group.by = 'integrated_clusters', label = T, label.size = 2) +
  theme(legend.position = 'none') +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5), byrow = TRUE)) + 
  ggtitle("Multi-species striatum datasets: integrated main cell types") &
  FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
  theme(legend.text=element_text(size=6), 
        strip.text.x = element_text(size = 5) , 
        legend.margin=margin(0,-20,0,0), legend.box.spacing = unit(0, "pt"),
        legend.box.margin=margin(0,-20,0,0),
        legend.spacing.x = unit(0, 'mm'), legend.spacing.y = unit(-2, 'mm')
  )


## make the main figure plot
pdf(here::here(PLOTDIR, 'multispecies_integrated_clusters_byProject.umap.pdf'),
    width = 7, height = 5)
p1 + p2 + plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))
dev.off()

## plot the subclusters, but only plot the gene expression 
cols4 = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                 c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                   'D2-Matrix', 'D2-Striosome', 'IC'))
obj@meta.data$integrated_subclusters = factor(obj$integrated_subclusters, names(cols4))
table(obj$integrated_subclusters)
p3 = DimPlot(obj, reduction = "umap", split.by = "Species", cols = cols4,
             group.by = 'integrated_subclusters', label = T, label.size = 2) +
  theme(legend.position = 'none') +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5), byrow = TRUE)) + 
  ggtitle("Multi-species striatum datasets: integrated subtypes") &
  FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
  theme(legend.text=element_text(size=6), 
        strip.text.x = element_text(size = 5) , 
        legend.margin=margin(0,-20,0,0), legend.box.spacing = unit(0, "pt"),
        legend.box.margin=margin(0,-20,0,0),
        legend.spacing.x = unit(0, 'mm'), legend.spacing.y = unit(-2, 'mm')
  )

pdf(here::here(PLOTDIR, 'multispecies_integrated_subclusters_bySpecies.umap.pdf'),
    width = 7, height = 3)
p3 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))
dev.off()


# Violin plot - Visualize single cell expression distributions in each cluster
obj2 = obj[,!obj$Project %in% c('Corces et al.', 'Li et al.')]
obj2@meta.data$Project = droplevels(obj2$Project)

features = c('DRD1', 'DRD2', 'DRD3', 'FOXP2')
pdf(here::here(PLOTDIR, 'gene_integrated_clusters.vln.pdf'), width = 7, height = 5)
VlnPlot(obj2, features = features, assay = 'SCT', group.by = 'integrated_clusters', 
        ncol = 1, pt.size = F, cols = cols3[-5]) & 
  facet_grid(~obj2$Project) & theme(axis.title.x = element_blank()) &
  coord_flip() & scale_x_discrete(limits = rev) &
  FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
  theme(strip.text.x = element_text(size = 5),
        axis.title.y = element_blank(), axis.text.x = element_text(angle = 0), 
        plot.title = element_blank()
        )
dev.off()


# Violin plot - Visualize single cell expression distributions in each cluster
pdf(here::here(PLOTDIR, 'gene_integrated_subclusters.vln.pdf'), width = 7, height = 6)
VlnPlot(obj2, features = features, assay = 'SCT', group.by = 'integrated_subclusters', 
        ncol = 1, pt.size = F, cols = cols4) & 
  facet_grid(~obj2$Project) & theme(axis.title.x = element_blank()) &
  coord_flip() & scale_x_discrete(limits = rev) &
  FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
  theme(strip.text.x = element_text(size = 5),
        axis.title.y = element_blank(), axis.text.x = element_text(angle = 0), 
        plot.title = element_blank()
  )
dev.off()


DefaultAssay(obj2) = 'RNA'
pdf(here::here(PLOTDIR, 'DRD_expression_byProject.umap.pdf'), width = 8.5, height = 10)
FeaturePlot(obj2, c("DRD1", "DRD2", 'DRD3', 'FOXP2'),  split.by = "Project") & 
  viridis::scale_color_viridis(option="magma")
dev.off()


pdf(here::here(PLOTDIR, 'QC_plots_byProject.umap.pdf'), width = 8.5, height = 5)
FeaturePlot(obj2, c("nCount_SCT", "nFeature_SCT"), split.by = "Project") & 
  viridis::scale_color_viridis(option="magma") &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()


pdf(here::here(PLOTDIR, 'QC_plots_byProject.vln.pdf'), width = 8, height = 4)
p1 = obj2[[]] %>% 
  ggplot(aes(y = nCount_RNA, x = integrated_clusters)) + 
  geom_violin(aes(fill = integrated_clusters)) + 
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = cols3, guide = 'none') +
  facet_grid( ~ Project) +  theme_bw(base_size = 8) + 
  ylab('# Detected UMI per cell')+ xlab('')

p2 = obj2[[]] %>% 
  ggplot(aes(y = nFeature_RNA, x = integrated_clusters)) + 
  geom_violin(aes(fill = integrated_clusters)) + 
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = cols3, guide = 'none') +
  facet_grid( ~ Project) + theme_bw(base_size = 8) +
  ylab('# Detected Genes per cell') + xlab('')

p1 + p2 + plot_layout(ncol = 1)+ 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))

dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_byProject.umap.pdf'), width = 8.5, height = 7)
DimPlot(obj, reduction = "umap", split.by = "Project", 
        cols = ArchR::paletteDiscrete(obj$Replicate), 
        group.by = 'Replicate') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()


# co-expression scatter plots
DefaultAssay(obj2) = 'SCT'

p1 = FeatureScatter(obj2, 'DRD1', 'DRD2', group.by = 'integrated_clusters', 
         cols = cols3, pt.size = 1) & 
  facet_grid(~obj2$Project, space = 'free') &
  theme(plot.title = element_blank(), legend.position = 'none', 
        strip.text.x = element_text(size = 8)) 

p2 = FeatureScatter(obj2, 'DRD1', 'DRD3',  group.by = 'integrated_clusters', 
               cols = cols3, pt.size = 1) & 
  facet_grid(~obj2$Project, space = 'free') &
  theme(plot.title = element_blank(), legend.position = 'none', 
        strip.text.x = element_text(size = 8)) 

p3 = FeatureScatter(obj2, 'DRD2', 'DRD3',  group.by = 'integrated_clusters', 
               cols = cols3, pt.size = 1) & 
  facet_grid(~obj2$Project, space = 'free') &
  theme(plot.title = element_blank(), legend.position = 'bottom', 
        strip.text.x = element_text(size = 8)) 

pdf(here::here(PLOTDIR, 'DRD_co-expression_byProject.scatter.pdf'), 
    width = 8, height = 8)

p1 + p2 + p3 + 
  plot_layout(ncol = 1)+ 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))
dev.off()




