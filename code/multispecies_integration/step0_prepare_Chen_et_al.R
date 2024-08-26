library(tidyverse)
library(Matrix)
library(data.table)
library(SeuratDisk)
library(Seurat)
library(here)

ss <- function(x, pattern, slot = 1, ...) {
  sapply(strsplit(x = x, split = pattern, ...), "[", slot)
}
DATADIR = here('data/tidy_data')
# export the object
obj = readRDS(here('data/tidy_data/Mouse_Chen/Chen2021_mouse_NAc_snRNA_msn.seurat.rds'))

# convert to human genes
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = LayerData(obj, layer = 'counts')
meta = obj[[]]
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts_hg)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]

obj2 = CreateSeuratObject(counts_hg, project = "Chen2021", assay = "RNA",
                          meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 

# save the object
save_fn2 = here(DATADIR,'Mouse_Chen/Chen2021_mouse_NAc_snRNA_msn.hgGenes.rds')
saveRDS(obj2, save_fn2)





