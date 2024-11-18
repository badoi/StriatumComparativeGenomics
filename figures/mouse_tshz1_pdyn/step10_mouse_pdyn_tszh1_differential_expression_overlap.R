library(tidyverse)
library(Seurat)
library(data.table)
library(patchwork)
library(cluster)
library(future)
library(here)
library(limma)
library(edgeR)
library(ggrepel)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
TABLDIR = 'figures/mouse_tshz1_pdyn/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(plot.title = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.5, "lines"))

##################################################################
## 1) read in the differential expression results for Pdyn and Tshz1
df_pdyn = here(TABLDIR, 'differential_expression_byPdyn_WMB-10Xv3-STR.NAc.xlsx') %>% 
  readxl::read_excel() %>% filter(Bonf.P.Val < 0.001, logFC > 0) %>% 
  mutate(group = 'Pdyn')

df_tshz1 = here(TABLDIR, 'differential_expression_byTshz1_WMB-10Xv3-STR.NAc.xlsx') %>% 
  readxl::read_excel() %>% filter(Bonf.P.Val < 0.001, logFC > 0) %>% 
  mutate(group = 'Tshz1')

df_deg_sig = bind_rows(df_pdyn, df_tshz1) %>% 
  group_by(genes) %>% mutate(tmp = n()) %>% ungroup() %>%
  mutate(diff_group = case_when(tmp > 1 ~ 'Shared\nMarkers', 
                                group == 'Tshz1' ~ 'Tshz1\nUnique',
                                group == 'Pdyn' ~ 'Pdyn\nUnique') %>% 
           factor(levels = c('Pdyn\nUnique',  'Tshz1\nUnique', 'Shared\nMarkers'))) %>% 
  group_by(group, diff_group) %>% tally() %>% arrange(diff_group) %>% 
  group_by(group) %>% mutate(label_y = cumsum(n) - n/2)

##################################################################
## 2) plot the integrated subclusters labels for the mouse NAcc
pdf(here(PLOTDIR, 'tshz1_pdyn_deg_overlap.barplot.pdf'), 
    height = 1.5, width = 1.5)
ggplot(df_deg_sig, aes(x = group, y = n, fill = diff_group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label = n, y = n), size = 2.5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#2D2D72','#b856d7', 'darkgray'), name = '') +
  labs(x = 'Correlated', y = 'Number of Marker Genes') +
  mytheme + theme(legend.position = 'bottom')
dev.off()

