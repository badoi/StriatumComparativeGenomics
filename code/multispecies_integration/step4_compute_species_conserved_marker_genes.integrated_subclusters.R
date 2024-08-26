library(tidyverse)
library(Seurat)
library(data.table)
library(limma)
library(edgeR)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

TABLDIR = 'figures/multispecies_integration/tables'
dir.create(TABLDIR, recursive = T)

source('https://github.com/YOU-k/voomByGroup/raw/main/voomByGroup.R')

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
## 1) create and load the pseudo bulk object for the analysis
pb_fn = here('data/tidy_data/rdas', 
             'integrated_multispecies_striatum_ARC_P15.pseudobulk.rds')
if(!file.exists(pb_fn)){
  obj = here('data/tidy_data/rdas', 
             'integrated_multispecies_striatum_ARC_P15.filtered.rds') %>% 
    readRDS()
  
  DefaultAssay(obj) = 'RNA'
  obj = JoinLayers(obj)
  
  ## grab the metadata for the pseudobulk object
  meta.data = obj[[]] %>% 
    group_by(Project, Replicate, Assay, Species, Sex, Case, Region, 
             integrated_clusters, integrated_subclusters) %>% 
    summarize(num_cells = n()) %>% 
    mutate(Replicate = str_replace_all(Replicate, '_', '-'),
           Region = str_replace_all(Region, '_', '-'),
           orig.ident = paste(Project, Region, Replicate, 
                              integrated_subclusters, sep = '_'))
  
  obj_pb = AggregateExpression(obj, assays = 'RNA', return.seurat = T,
                               group.by = c('Project', 'Region', 'Replicate', 
                                            'integrated_subclusters'))
  
  meta.data = meta.data[match(obj_pb$orig.ident, meta.data$orig.ident),]
  obj_pb = AddMetaData(obj_pb, meta.data)
  
  ## save the pseudobulk object
  obj_pb %>% saveRDS()
} else {
  obj_pb = readRDS(pb_fn)
}


###########################################
## 2) limma-voom cell type marker genes 
obj_pb = obj_pb %>% subset(Assay == 'RNA')
obj_pb$cdr = (LayerData(obj_pb, 'counts') > 0) %>% colMeans() %>% scale()

y <- DGEList(counts = LayerData(obj_pb, 'counts'), genes = rownames(obj_pb))
A <- rowSums(y$counts)
isexpr <- A > 50

hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) # 36079  4518

# Apply scale normalization:
y <- calcNormFactors(y)


# Create a design matrix for this experiment as a hypothesis of the gene expression
meta2 =  obj_pb[[]] %>% mutate(Project = make.names(Project), 
         integrated_subclusters = make.names(integrated_subclusters))
design <- model.matrix(~integrated_subclusters + Project + Sex + Region + Case + 
                         # the last are the known covariates with pseudobulk data
                          num_cells + cdr + 0,  data = meta2)

# extract the top differentially expressed genes, for example for 
design2 = design
colnames(design2) = colnames(design) %>% make.names()

celltypes = colnames(design2)[1:7] %>% str_remove('integrated_subclusters') %>%
  set_names()

## make the cell type contrasts in Putamens
con_celltypes = sapply(celltypes, function(cell) {
  fgd = colnames(design2)[1:7] %>% str_subset(cell) 
  bgd = colnames(design2)[1:7] %>% str_subset(cell, negate = T) 
  N_bgd = bgd %>% length()
  bgd = bgd %>% paste(collapse = ' + ')
  
  ## export the contrast
  paste('(',fgd,') - (',bgd,')/',N_bgd)
})

cont.matrix <- makeContrasts(contrasts= con_celltypes, levels=design2)
rownames(cont.matrix) = colnames(design)

# Fit the voom model, using the sophisticated quality weights to account for
# sample-wise variability in the pseudobulk data, 
v <- voomByGroup(y, design=design, group= obj_pb$Replicate)
vfit <- lmFit(v) %>% eBayes()

fit2 <- contrasts.fit(vfit, cont.matrix) %>% eBayes()

## compute the DEGs from these contrasts
df_degs = lapply(setNames(colnames(cont.matrix), celltypes), topTable, 
                  fit = fit2, n = Inf) %>% rbindlist(idcol = 'celltype') %>% 
  mutate(adj.P.Val = p.adjust(P.Value, method = 'fdr'), 
         Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni')) %>% 
  filter(Bonf.P.Val < 0.01, logFC > 1)

df_degs %>% count(celltype) 

#        celltype     n
# 1:       D1.D2H   250
# 2:    D1.Matrix   116
# 3:     D1.NUDAP  1204
# 4: D1.Striosome    70
# 5:    D2.Matrix   108
# 6: D2.Striosome    79
# 7:           IC   434

df_degs %>% filter( grepl('PDYN|OPRM1', genes))
df_degs %>% filter(genes == 'TSHZ1')
df_degs %>% filter(grepl('^DRD', genes))
df_degs %>% filter(genes == 'DRD3')

## write out the DEGs
out_fn = here(TABLDIR, 'integrated_multispecies_striatum.integrated_subclusters.xlsx')
df_degs %>% writexl::write_xlsx(out_fn)


