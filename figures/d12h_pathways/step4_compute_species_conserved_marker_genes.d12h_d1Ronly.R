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
PLOTDIR = 'figures/d12h_pathways/plots'
dir.create(PLOTDIR, recursive = T)

TABLDIR = 'figures/d12h_pathways/tables'
dir.create(TABLDIR, recursive = T)

source('https://github.com/YOU-k/voomByGroup/raw/main/voomByGroup.R')
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#############################################################
## 1) create and load the pseudo bulk object for the analysis
pb_fn = here('data/tidy_data/rdas', 
             'integrated_multispecies_striatum_ARC_P15.pseudobulk.rds')

obj_pb = readRDS(pb_fn)
obj_pb = obj_pb %>% subset(integrated_clusters %in% c('D1') | 
                             integrated_subclusters %in% c('D1/D2H')) %>%
  subset(Assay == 'RNA')
obj_pb$cdr = (LayerData(obj_pb, 'counts') > 0) %>% colMeans() %>% scale()

table(obj_pb$integrated_clusters, obj_pb$integrated_subclusters)

###########################################
## 2) limma-voom cell type marker genes 
y <- DGEList(counts = LayerData(obj_pb, 'counts'), genes = rownames(obj_pb))
A <- rowSums(y$counts)
isexpr <- A > 50

hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) #  34064  2031

# Apply scale normalization:
y <- calcNormFactors(y)

# Create a design D1.2R only for this experiment as a hypothesis of the gene expression
design <- model.matrix(~integrated_clusters + Project + Sex + Case + 
                         # the last are the known covariates with pseudobulk data
                         num_cells + cdr,  data = obj_pb[[]])

# fit the linear model
v <- voomByGroup(y, design=design, group= obj_pb$Replicate)
vfit <- lmFit(v) %>% eBayes()

## compute the DEGs from these contrasts
df_degs = topTable(fit = vfit, coef = 'integrated_clusterseSPN', n = Inf) %>% 
  mutate(Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni'), 
         Marker = ifelse(logFC > 0, 'D12H', 'D1R-only')) %>% 
  rename_at(vars(logFC, t), ~paste0(.x, '_D12H_vs_D1R-only'))

## write out the DEGs
out_fn = here(TABLDIR, 'integrated_multispecies_striatum.D1R-onlyD12H.xlsx')
df_degs %>% writexl::write_xlsx(out_fn)


#####################################################
## 3) limma-voom cell type marker genes in rodents
obj_pb_mm = obj_pb %>% subset(Species %in% c('Mouse', 'Rat') )
y <- DGEList(counts = LayerData(obj_pb_mm, 'counts'), genes = rownames(obj_pb_mm))
A <- rowSums(y$counts)
isexpr <- A > 50
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) #  14943  1050

# Apply scale normalization:
y <- calcNormFactors(y)
design <- model.matrix(~integrated_clusters + Project + Sex + Case + 
                         num_cells + cdr,  data = obj_pb_mm[[]])
vfit <- voomByGroup(y, design=design, group= obj_pb_mm$Replicate) %>% 
  lmFit() %>% eBayes()

## compute the DEGs from these contrasts
df_deg_mm = topTable(fit = vfit, coef = 'integrated_clusterseSPN', n = Inf) %>% 
  mutate(Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni'), 
         Marker = ifelse(logFC > 0, 'D12H', 'D1R-only')) %>% 
  rename_at(vars(logFC, t), ~paste0(.x, '_D12H_vs_D1R-only')) 

## write out the DEGs
out_fn = here(TABLDIR, 'integrated_multispecies_striatum.D1R-onlyD12H.rodent.xlsx')
df_deg_mm %>% writexl::write_xlsx(out_fn)



#####################################################
## 4) limma-voom cell type marker genes in primates
obj_pb_hg = obj_pb %>% subset(Species %in% c('Human', 'Macaque') )
y <- DGEList(counts = LayerData(obj_pb_hg, 'counts'), genes = rownames(obj_pb_hg))
A <- rowSums(y$counts)
isexpr <- A > 50
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) #  14943  1050

# Apply scale normalization:
y <- calcNormFactors(y)
design <- model.matrix(~integrated_clusters + Project + Sex + Case + 
                         num_cells + cdr,  data = obj_pb_hg[[]])
vfit <- voomByGroup(y, design=design, group= obj_pb_hg$Replicate) %>% 
  lmFit() %>% eBayes()

## compute the DEGs from these contrasts
df_deg_hg = topTable(fit = vfit, coef = 'integrated_clusterseSPN', n = Inf) %>% 
  mutate(Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni'), 
         Marker = ifelse(logFC > 0, 'D12H', 'D1R-only')) %>% 
  rename_at(vars(logFC, t), ~paste0(.x, '_D12H_vs_D1R-only')) 

## write out the DEGs
out_fn = here(TABLDIR, 'integrated_multispecies_striatum.D1R-onlyD12H.primate.xlsx')
df_deg_hg %>% writexl::write_xlsx(out_fn)

