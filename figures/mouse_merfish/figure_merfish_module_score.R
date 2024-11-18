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
PLOTDIR = 'figures/mouse_merfish/plots'
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
hg38_to_mm10 = here('/projects/pfenninggroup/machineLearningForComputationalBiology',
                    'resources/genomes/GRCh38.p13',
                    'ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv') %>% 
  read_tsv() %>% rename_with(make.names) %>% 
  filter(!duplicated(Mouse.gene.name), !is.na(Mouse.gene.name)) %>% 
  select(Gene.name, Mouse.gene.name) %>% deframe()

markers_dv = here('figures/multispecies_integration/tables',
                  'integrated_multispecies_striatum.dorsalVentral.xlsx') %>% 
  readxl::read_xlsx() %>% rename('Marker' = 'RegionMarker')

markers_sm = here('figures/multispecies_integration/tables',
                  'integrated_multispecies_striatum.MatrixStriosome.xlsx') %>% 
  readxl::read_xlsx() %>% rename('Marker' = 'CompartmentMarker')

markers_df = bind_rows(markers_dv, markers_sm) %>% 
  mutate(mouse_genes = hg38_to_mm10[genes]) %>% filter(!is.na(mouse_genes)) %>% 
  filter(Bonf.P.Val < 0.05) 

regions = c("Dorsal" , 'Ventral', 'Striosome', 'Matrix')

##############################################################################
## 2) load the merfish datasets and compute the module score and plot
str_merf_fn = here(DATADIR) %>%
  list.files(pattern = 'integrated_striatum.seurat.rds', full.names = T)

i = 1
for (i in seq_along(str_merf_fn)) {
  obj = readRDS(str_merf_fn[i])
  
  if (! all(regions %in% names(obj[[]])) ) {
    ## compute the module score for each cluster
    markers_list = markers_df %>% filter(mouse_genes %in% Features(obj)) %>%
      split(x = .$mouse_genes, f = .$Marker)
    
    cells_rankings <- AUCell_buildRankings(LayerData(obj, layer = 'data'), plotStats=F)
    cells_AUC <- AUCell_calcAUC(markers_list, cells_rankings, nCores = 4)
    df_module_score = cells_AUC %>% getAUC() %>% t() %>% as.data.frame()
    df_module_score = df_module_score[match(Cells(obj), rownames(df_module_score)),]
    
    obj = obj %>% AddMetaData(df_module_score)
    obj %>% saveRDS(str_merf_fn[i])
  }
}



##############################################################################
## 3) plot the module score for Striosome-Matrix and Ventral-Dorsal on MERFISH
for (i in seq_along(str_merf_fn) %>% rev()) {
  obj = readRDS(str_merf_fn[i])
    
  ## find the images that have both ACB and CP regions to plot the module scores
  table(obj$parcellation_substructure)
  to_plot_images = obj[[]] %>% group_by(brain_section_label) %>% 
    filter(all(c('ACB','CP') %in% parcellation_substructure)) %>% ungroup() %>% 
    pull(brain_section_label) %>% unique()
  
  ## make plots to evaluate the transferred integrated cluster labels
  for(fov in to_plot_images){
    plot_fn = str_merf_fn[i] %>% basename() %>% 
      str_remove('integrated_striatum.seurat.rds') %>% 
      paste0(., fov, '.regionModule.pdf') %>% here(PLOTDIR,'regionModule',.)
    dir.create(dirname(plot_fn), recursive = T, showWarnings = F)
    # if ( file.exists(plot_fn) ) next
    pdf(plot_fn, height = 3, width = 1)
    pp = obj %>% 
      SpatialFeaturePlot(features = regions, pt.size.factor = 1.5, 
                         alpha = c(0.2, 1), ncol = 1, max.cutoff = 'q95', 
                         min.cutoff = 'q5', images = make.names(fov)) & 
      # DarkTheme() & 
      NoGrid() & NoAxes() & theme(legend.position = "none")
    print(pp)
    dev.off()
  }
}


fov = 'C57BL6J-638850.52'
plot_fn = str_merf_fn[i] %>% basename() %>% 
  str_remove('integrated_striatum.seurat.rds') %>% 
  paste0(., '.legend.pdf') %>% here(PLOTDIR,'regionModule',.)

pdf(plot_fn, height = .75, width = 4)
pp = obj %>% 
  SpatialFeaturePlot(features = regions, pt.size.factor = 1.5, 
                     alpha = c(0.2, 1), ncol = 4, max.cutoff = 'q95', 
                     min.cutoff = 'q5', images = make.names(fov)) & 
  # DarkTheme() & 
  guides(fill = guide_colorbar(barwidth = .25, barheight = 2)) &
  NoGrid() & NoAxes() & theme(legend.position = "right", ,
                              legend.title=element_text(size=4), 
                              legend.text=element_text(size=4), 
                              legend.key.width  = unit(0.5, "lines"),
                              legend.key.height = unit(1, "lines"))
print(pp)
dev.off()


