library(Seurat)
library(tidyverse)
library(here)

# Load the Seurat object to split
source_fn = here('data/tidy_data/Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.rds')
obj_mouse = readRDS(source_fn)

