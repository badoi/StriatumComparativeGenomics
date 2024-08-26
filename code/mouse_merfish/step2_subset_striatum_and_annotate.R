library(tidyverse)
library(here)
library(Seurat)

DATADIR='data/tidy_data/mouse_merfish'
dir.create(here(DATADIR,'plots'), showWarnings = F)

############################
## 0) prepare the metadata
clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)

subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

#####################################################################
## 1) Load the data and convert Allen cluster name to cell type name
obj = readRDS(here('data/tidy_data/rdas/integrated_rodent_striatum.supercells.rds'))

obj = obj %>% subset(Project == 'Zeng et al.') %>% 
  FindVariableFeatures(verbose = F) %>% RunPCA(verbose = F) %>% 
  RunUMAP(dims = 1:30, return.model = T, verbose = F)

##############################################################################
## 2) read in each MERFISH dataset and transfer the integrated clusters label
orig_merf_fn = here(DATADIR) %>%
  list.files(pattern = 'spatial.seurat.rds', full.names = T)
str_merf_fn = orig_merf_fn %>% 
  str_replace('spatial.seurat.rds', 'integrated_striatum.seurat.rds')

for ( i in seq_along(orig_merf_fn) ) {
  if ( !file.exists(str_merf_fn[i]) ) {
    obj_merf = readRDS(orig_merf_fn[i]) 
    
    class_keep = grep('CNU GABA|CNU-LGE GABA', obj_merf$class, value = T) %>% unique()

    ## subset the data to the striatum and MSNs
    obj_merf = obj_merf %>% subset(parcellation_division == 'STR') %>% 
      subset(class %in% class_keep) %>% ScaleData(verbose = F) %>% 
      FindVariableFeatures(verbose = F) %>% RunPCA(assay = "RNA", verbose = F) 
    
    ## map the integrated cross-species striatum subtype to the MERFISH data
    anchors <- FindTransferAnchors(reference = obj, query = obj_merf, k.anchor = 25, 
                                   normalization.method = "LogNormalize", )
    
    integrated_clusters <- 
      TransferData(anchors, refdata = obj$integrated_clusters,
                   weight.reduction = obj_merf[["pca"]], dims = 1:30)
    
    integrated_subclusters <-
      TransferData(anchors, refdata = obj$integrated_subclusters,
                   weight.reduction = obj_merf[["pca"]], dims = 1:30)
    
    obj_merf$integrated_clusters <- integrated_clusters$predicted.id
    obj_merf$integrated_subclusters <- integrated_subclusters$predicted.id
    
    ## save the transferred data
    obj_merf %>% saveRDS(str_merf_fn[i])
  }
}


#########################################################################
## 3) plot the spatial arrangement of the integrated labels MERFISH data
for(file in str_merf_fn){
  obj_merf = readRDS(file)
  
  ## make plots to evaluate the transferred integrated cluster labels
  for(fov in unique(obj_merf$brain_section_label)){
    plot_fn = here(DATADIR,'plots', paste0(fov, ".integrated_subclusters.pdf"))
    if ( file.exists(plot_fn) ) next
    ## plot the spatial data
    pdf(plot_fn)
    pp = obj_merf %>% 
      SpatialDimPlot(group.by = "integrated_subclusters", pt.size.factor = 2,
                     cols = cols2, images = make.names(fov)) & 
      theme(legend.position = "bottom")
    print(pp)
    dev.off()
  }
}


