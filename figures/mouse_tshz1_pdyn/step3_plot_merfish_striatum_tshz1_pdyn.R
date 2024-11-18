library(tidyverse)
library(AUCell)
library(Seurat)
library(data.table)
library(patchwork)
library(viridis)
library(cluster)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/mouse_merfish'
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
dir.create(PLOTDIR, recursive = T)

#############################
## 0) prepare the metadata
clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)
subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

##############################################################################
## 1) read in ventral-dorsal, striosome-matrix markers to compute module score
str_merf_fn = here(DATADIR) %>%
  list.files(pattern = 'integrated_striatum.seurat.rds', full.names = T)

for (i in seq_along(str_merf_fn) %>% rev()) {
  obj = readRDS(str_merf_fn[i])
  
  ## find the images that have both ACB and CP regions to plot the module scores
  to_plot_images = obj[[]] %>% group_by(brain_section_label) %>% 
    filter(all(c('ACB') %in% parcellation_substructure)) %>% ungroup() %>% 
    pull(brain_section_label) %>% unique()
  
  obj_list = list('eSPN' =  obj %>% subset(integrated_clusters == 'eSPN'), 
                  'D1' = obj %>% subset(integrated_clusters == 'D1'))
  
  ## make plots to evaluate the transferred integrated cluster labels
  height = 1.25
  genes = c('Tshz1', 'Pdyn', 'Oprm1')
  genes_plot = genes[genes %in% Features(obj)]
  
  for(fov in to_plot_images) {
    ## plot the integrated clusters on the MERFISH data
    plot_fn = str_merf_fn[i] %>% basename() %>% 
      str_remove('integrated_striatum.seurat.rds') %>% 
      paste0(., fov, '.integrated_subclusters.pdf') %>% 
      here(PLOTDIR, 'merfish_integrated_clusters',.)
    dir.create(dirname(plot_fn), recursive = T, showWarnings = F)
    pdf(plot_fn, height = height, width = 1)
    pp = obj %>% 
      SpatialDimPlot(pt.size.factor = 3, images = make.names(fov), 
                     group.by = 'integrated_subclusters', cols = cols2) & 
      NoGrid() & NoAxes() & theme(legend.position = "none")
    print(pp)
    dev.off()
    
    ## plot the expression of target genes measured by the MERFISH data
    plot_fn = str_merf_fn[i] %>% basename() %>% 
      str_remove('integrated_striatum.seurat.rds') %>% 
      paste0(., fov,'.',paste(genes_plot, collapse = '-'),'.pdf') %>%
      here(PLOTDIR, 'merfish_expression',.)
    dir.create(dirname(plot_fn), recursive = T, showWarnings = F)
    
    pdf(plot_fn, height = height * length(genes_plot), width = 1)
    pp = obj %>% 
      SpatialFeaturePlot(features = genes_plot, , pt.size.factor = 4, 
                         alpha = c(0.2, 1), ncol = 1, max.cutoff = 'q95', 
                         min.cutoff = 'q5', images = make.names(fov)) & 
      NoGrid() & NoAxes() & theme(legend.position = "bottom", ,
                                  legend.title = element_text(size=4), 
                                  legend.text = element_text(size=4), 
                                  legend.key.width  = unit(0.5, "lines"),
                                  legend.key.height = unit(0.25, "lines"))
    
    print(pp)
    dev.off()
  }
}

