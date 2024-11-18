library(tidyverse)
library(treemap)
library(readxl)
library(here)
library(ggh4x)

PLOTDIR='figures/d12h_pathways/plots'
DATADIR='figures/d12h_pathways/tables'
dir.create(here(PLOTDIR), recursive = T)
dir.create(here(DATADIR), recursive = T)

#
mapping_df = 
  read_xlsx(here(DATADIR, 'syngo_ontologies_with_annotations_matching_user_input.xlsx')) %>% 
  rename_with(make.names)


d1h_df = here(DATADIR, 'SynGO_geneset_analysis__D1-2H__2024-07-06 17;54', 'syngo_ontologies_with_annotations_matching_user_input.xlsx') %>% 
 read_xlsx() %>% rename_with(make.names) %>% 
  filter(!is.na(GSEA..gene.cluster..FDR.corrected.p.value)) %>%
  filter(!is.na(genes...without.hierarchical.rollup...hgnc_symbol)) %>%
  mutate(celltype = 'D1/D2H')

d1_df = here(DATADIR, 'SynGO_geneset_analysis__D1R-only__2024-07-06 17;57', 'syngo_ontologies_with_annotations_matching_user_input.xlsx') %>%
  read_xlsx() %>% rename_with(make.names) %>% 
  filter(!is.na(GSEA..gene.cluster..FDR.corrected.p.value)) %>%
  filter(!is.na(genes...without.hierarchical.rollup...hgnc_symbol)) %>%
  mutate(celltype = 'D1R-only')

group_label = d1_df %>% 
  mutate(group = str_sub(user.interface.reference.code, 1, 1) %>% paste0(GO.domain),
         group_number = str_sub(user.interface.reference.code, 2,2)) %>% 
  filter(group_number == '1') %>% 
  select(group, GO.term.name) %>% deframe()

enrich_df = bind_rows(d1h_df, d1_df) %>% 
  mutate(celltype = factor(celltype, levels = c('D1/D2H', 'D1R-only'))) %>% 
  filter(!is.na(GO.parent.term.ID)) %>% 
  mutate(label = genes...without.hierarchical.rollup...hgnc_symbol %>% 
           str_split(';') %>% map_chr(~sort(.x) %>% paste(collapse = ', ')),
         GO.term.name = GO.term.name %>% 
           str_replace_all('neurotransmitter', 'NT-') %>% 
           str_replace_all('transmitter', 'NT') %>% 
           str_replace_all('receptor', 'R') %>%
           str_replace_all('membrane', 'mem.') %>%
           str_replace_all('channel', 'chan.') %>%
           str_replace_all('involved in regulation of', 'regulating') %>%
           str_replace_all('synaptic', 'syn.'),
         values = ifelse(is.na(GSEA..gene.cluster..FDR.corrected.p.value), 0,
                         -log10(GSEA..gene.cluster..FDR.corrected.p.value)), 
         values = ifelse(celltype == 'D1/D2H', -values, values), 
         group = str_sub(user.interface.reference.code, 1, 1) %>% paste0(GO.domain), 
         group_num = str_sub(user.interface.reference.code, 2) %>% as.numeric(), 
         group_name = group_label[group] %>% 
           str_replace_all('neurotransmitter', 'NT-') %>% 
           str_replace_all('transmitter', 'NT') %>% 
           str_replace_all('receptor', 'R') %>%
           str_replace_all('membrane', 'mem.') %>%
           str_replace_all('channel', 'chan.') %>%
           str_replace_all('involved in regulation of', 'regulating') %>%
           str_replace_all('synaptic', 'syn.'), 
         yintercept =ifelse(celltype == 'D1/D2H', log10(0.05), -log10(0.05)),
         hjust =ifelse(celltype == 'D1/D2H', 1, 0)) %>% 
  filter(group_num < 26) %>% 
  group_by(GO.term.name) %>% mutate(tmp = mean(abs(values))) %>%
  ungroup() %>% arrange(group_name, tmp) %>% 
  mutate(GO.term.name = factor(GO.term.name, levels = unique(GO.term.name)))

p1 = enrich_df %>% 
  ggplot(aes(x = GO.term.name, y = values, fill = celltype)) + 
  geom_hline(aes(yintercept = yintercept), linetype = 2, color = 'red', 
             size =.25) +
  geom_col() + 
  geom_text(aes(label = label, hjust = hjust), y = 0, vjust = 0.5, size = .6) +
  facet_grid2(GO.domain + group_name ~ celltype, scales = 'free', 
             space = 'free_y', axes = 'y', switch = 'y') + 
  coord_flip() + 
  scale_y_continuous(labels = abs) +
  theme_bw(base_size = 3) + 
  labs(y = '-log10(FDR)', x = 'Synaptic Onotology') +
  scale_fill_manual(values = c('D1/D2H' = '#33a02c', 'D1R-only' = '#1f78b4')) +
  theme(axis.text.y = element_text(hjust = 0.5), legend.position = 'none', 
        strip.background = element_rect(fill = NA, color = "black"), 
        strip.placement = "outside", 
        panel.spacing.y = unit(.25, "lines")) 
p1

pdf(here(PLOTDIR, 'syngo_enrichment_d1h_d1r.pdf'), width = 4, height = 2)
p1
dev.off()


