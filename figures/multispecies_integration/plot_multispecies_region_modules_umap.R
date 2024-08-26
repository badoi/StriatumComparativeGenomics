library(tidyverse)
library(Seurat)
library(data.table)
library(patchwork)
library(cluster)
library(future)
library(here)
library(AUCell)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

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
multispecies_fn = 
  here('data/tidy_data/rdas/integrated_multispecies_striatum_ARC_P15.filtered.rds')
obj = readRDS(multispecies_fn)

DefaultAssay(obj) = 'RNA'

if (!'counts' %in% Layers(obj)) {
  obj = JoinLayers(obj)
  obj %>% saveRDS(multispecies_fn)
}

obj = obj %>% FindVariableFeatures() %>% ScaleData()

### read in the ventral-dorsal, striosome-matrix markers to compute module score
regions = c("Dorsal" , 'Ventral', 'Striosome', 'Matrix')
markers_dv = here('figures/multispecies_integration/tables',
                  'integrated_multispecies_striatum.dorsalVentral.xlsx') %>% 
  readxl::read_xlsx() %>% rename('Marker' = 'RegionMarker')

markers_sm = here('figures/multispecies_integration/tables',
                  'integrated_multispecies_striatum.MatrixStriosome.xlsx') %>% 
  readxl::read_xlsx() %>% rename('Marker' = 'CompartmentMarker')

markers_df = bind_rows(markers_dv, markers_sm) %>% filter(Bonf.P.Val < 0.05) 
marker_list = markers_df %>% split(genes, f = .$Marker)
marker_list = marker_list[regions]

table(markers_df %in% Features(obj))

## compute the module score for each cluster
obj = AddModuleScore(obj, marker_list, assay = 'RNA', name = names(marker_list))

#############################################################
## 3) plot the module scores
pdf(here(PLOTDIR, 'region_module_scores.umap.pdf'), height = 2.75, width = 1)
FeaturePlot(obj, features = paste0(regions, 1:4), raster = T,  
            reduction = "umap.full", ncol = 1, cols = c('#FFFFFF00', '#000000FF'),
            min.cutoff = 'q5', max.cutoff = 'q95') & 
  FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  guides(fill = guide_colorbar(barwidth = .15, barheight = .7)) &
  theme(legend.position = 'right', title = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.line = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"), 
        )
dev.off()



#################################
## 4) plot the dopamine receptors
pdf(here(PLOTDIR, 'dopamine_receptors_feature.umap.pdf'), height = .75, width = 2.75)
FeaturePlot(obj, features = c('DRD1', 'DRD2', 'DRD3', 'FOXP2'), raster = T,  
            reduction = "umap.full", ncol = 4, cols = c('#FFFFFF', '#000000'),
            min.cutoff = 'q5', max.cutoff = 'q50') & 
  FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  guides(fill = guide_colorbar(barwidth = 0.05, barheight = .7)) &
  theme(legend.position = 'right', title = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.line = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"), 
  )
dev.off()







