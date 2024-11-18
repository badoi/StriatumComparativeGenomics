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
## 1) load the single cell embedding across multiple species
supercells_fn = 
  here('data/tidy_data/Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.supercells.rds')
obj_super = readRDS(supercells_fn) %>% subset(Region %in% 'NAcc') 

source('https://github.com/YOU-k/voomByGroup/raw/main/voomByGroup.R')

## extract the gene expression matrix for limma analysis
y <- DGEList(counts = LayerData(obj_super, 'counts'), genes = rownames(obj_super))
A <- rowSums(y$counts)
isexpr <- A > 5
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y) # 19317  5289

# Apply scale normalization:
y <- calcNormFactors(y)

obj_super$Pdyn = getCounts(y)['Pdyn',]

design <- model.matrix(~ Pdyn + Sex + size,  data = obj_super[[]])
vfit <- voomByGroup(y, design=design, group= obj_super$Replicate) %>% 
  lmFit() %>% eBayes()

## compute the DEGs from these contrasts
striosome_markers = c('EphaA5', 'Htr2a', 'Rxrg', 'Cdh11', 'Lgi1', 'Col25a1', 
                      'Phactr1', 'Myo3b', 'Ripor2', 'Khdrbs2', 'Nnat', 'Lypd1', 
                      'Pnoc', 'Nts', 'Kremen1', 'Tac1', 'Tacr1', 'Chat')
df_deg = topTable(fit = vfit, coef = 'Pdyn', n = Inf) %>% 
  mutate(Bonf.P.Val = p.adjust(P.Value, method = 'bonferroni'), 
         to_label = ifelse(genes %in% c('Oprm1', 'Pdyn', 'Tshz1', 
                                        striosome_markers), genes, NA), 
         diff_color = ifelse(genes %in% c('Oprm1', 'Pdyn', 'Tshz1'), T, F))

## write out the DEGs
out_fn = here(TABLDIR, 'differential_expression_byPdyn_WMB-10Xv3-STR.NAc.xlsx')
df_deg %>% writexl::write_xlsx(out_fn)


##################################################################
## 2) plot the integrated subclusters labels for the mouse NAcc
mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', plot.title = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

df_deg_sig = df_deg %>% filter(Bonf.P.Val < 0.001, logFC > 0)

pdf(here(PLOTDIR, 'pdyn_differential_expression.volcano.pdf'), 
    height = 1.5, width = 1.5)
ggplot(df_deg_sig, aes(x = logFC, y = -log10(Bonf.P.Val))) + 
  geom_point(aes(shape = logFC > 1, alpha = !is.na(to_label), color = diff_color),
             size = 1) +
  geom_text(aes(label = to_label), vjust = 'inward', 
            hjust = 'inward', size = 1.8) +
  scale_color_manual(values = c('#6d6dc7', '#2D2D72')) +
  scale_shape_manual(values = c(20, 19)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  labs(x = 'Log2(Correlation with Pdyn)', 
       y = '-log10(Bonferroni P-value)') +
  mytheme
dev.off()

