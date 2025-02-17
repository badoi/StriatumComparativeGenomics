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


#############################################################
## 1) create and load the pseudo bulk object for the analysis
pb_fn = here('data/tidy_data/rdas', 
             'integrated_multispecies_striatum_ARC_P15.pseudobulk.rds')

obj_pb = readRDS(pb_fn)
obj_pb = obj_pb %>% subset(Assay == 'RNA')
obj_pb$cdr = (LayerData(obj_pb, 'counts') > 0) %>% colMeans() %>% scale()
obj_pb$Region2 = ifelse(obj_pb$Region == 'NAcc', 'Ventral','Dorsal')


###########################################
## 2) limma-voom cell type marker genes 
y <- DGEList(counts = LayerData(obj_pb, 'counts'), genes = rownames(obj_pb))
A <- rowSums(y$counts)
isexpr <- A > 50

hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) # 36079  4518

# Apply scale normalization:
y <- calcNormFactors(y)

# Create a design matrix for this experiment as a hypothesis of the gene expression
design <- model.matrix(~Region2 + integrated_subclusters + Project + Sex + Case + 
                         # the last are the known covariates with pseudobulk data
                          num_cells + cdr,  data = obj_pb[[]])

# fit the linear model
v <- voomByGroup(y, design=design, group= obj_pb$Replicate)
vfit <- lmFit(v) %>% eBayes()

## compute the DEGs from these contrasts
df_degs = topTable(fit = vfit, coef = 'Region2Ventral', n = Inf) %>% 
  mutate(Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni'), 
         RegionMarker = ifelse(logFC > 0, 'Ventral', 'Dorsal')) %>% 
  rename_at(vars(logFC, t), ~paste0(.x, '_Ventral_vs_Dorsal')) %>% 
  filter(Bonf.P.Val < 0.05)

df_degs %>% count(RegionMarker) 

#   RegionMarker   n
# 1       Dorsal 138
# 2      Ventral 201

## write out the DEGs
out_fn = here(TABLDIR, 'integrated_multispecies_striatum.dorsalVentral.xlsx')
df_degs %>% writexl::write_xlsx(out_fn)


