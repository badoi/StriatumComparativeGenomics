library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
TABLDIR = 'figures/mouse_tshz1_pdyn/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

save_fn = here(DATADIR, 'Mouse_Zeng/WMB-10Xv3-STR.neuron.seurat.rds')

if (!file.exists(save_fn)) {
  obj = readRDS(here(DATADIR, 'Mouse_Zeng/WMB-10Xv3-STR.seurat.rds'))
  
  head(obj,2)
  
  table(obj$division)
  table(obj$region_of_interest_acronym)
  table(obj$division, obj$region_of_interest_acronym)
  table(obj$class, obj$region_of_interest_acronym)
  
  ind_neur = which(obj$region_of_interest_acronym %in% c('STRd', 'STRv') & 
                     !obj$class %in% c('07 LSX GABA', '12 MOB-CR Glut',
                                       '01 IT-ET Glut', '02 NP-CT-L6b Glut') &
                     obj$division %in% c('2 Subpallium GABAergic', 
                                         '4 CBX-MOB-other neuronal'))
  
  obj_neur = obj[,ind_neur] 
  
  ## subset to just MSNs from unaffected individuals
  saveRDS(obj_neur, save_fn)
} else {
  obj_neur = readRDS(save_fn)
}

obj_neur$subclass2 = str_remove(obj_neur$subclass, "^\\d+\\s*")
obj_neur$group = case_when(obj_neur$subclass %in% 
               c('051 MSN D1 Gaba', '052 MSN D2 Gaba', 
                 '050 OT D3 Folh1 Gaba', '053 MSN D1 Sema5a Gaba', 
                 '054 STR-PAL Chst9 Gaba') ~ 'MSNs',
             TRUE ~ "probably interneurons")

table(obj_neur$subclass2) %>% sort()
table(obj_neur$subclass2, obj_neur$group) 

# drop the subclasses with less than 10 cells
keep = names(table(obj_neur$subclass2))[table(obj_neur$subclass2) > 10]
obj_neur = obj_neur[,obj_neur$subclass2 %in% keep]

# Plotting
mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', 
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

genes_hs = c("Oprm1", "Pdyn", "Tshz1" ,"Tac1" ,"Tacr1", "Chat" ,
             "Nts", "Kremen1" ,"Pnoc")
cols = c(DiscretePalette(15, palette = 'alphabet', shuffle = T), 
         DiscretePalette(16, palette = 'alphabet2', shuffle = T))

pdf(here(PLOTDIR, 'allen_mouse_striatal_neurons_pdyn_tshz1_oprm1.vln.pdf'), 
    height = 6, width = 6)
VlnPlot(obj_neur, genes_hs, group.by = 'subclass2', cols = cols, 
        raster = T, pt.size = 0, same.y.lims = T) & mytheme & 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) & 
  facet_grid(~obj_neur$group, scales = 'free_x', space = 'free_x')
dev.off()

