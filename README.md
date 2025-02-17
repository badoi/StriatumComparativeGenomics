# **Multi-dimensional characterization of striatal projection neuron heterogeneity in the adult striatum**
by _Jenesis Gayden, Stephanie Puig, Chaitanya Srinivasan, BaDoi N. Phan, Ghada Abdelhady, Silas A. Buck, Mackenzie C. Gamble, Hugo A. Tejeda, Yan Dong, Andreas R. Pfenning, Ryan W. Logan, Zachary Freyberg_*

*Corresponding author:\
Zachary Freyberg, M.D., Ph.D.\
3811 O’Hara Street, BST W1640\
Pittsburgh, PA 15213\
Tel: (412) 648-0033 \
Email: freyberg at pitt dot edu

[Link to paper](https://www.biorxiv.org/content/10.1101/2023.05.04.539488v1)

## Multi-species integration of single cell striatal projection neurons
This repository includes the code for comparative integrative analyses of single cell datasets of the mammalian striatum. 

The main R-scripts to perform the integration and marker gene analyses are include in the `code/` directory.
- `step0_*.R`: performs dataset specific manipulations to subset data if not already down to striatal projection neurons based on original publication labels
- `step1_integrate_multispecies_striatum.ARC.R`: performs Seurat reciprocal PCA integration of the SCTransform normalized counts from snRNA-seq and snATAC-seq data in a guided manner
- `step2_compute_species_specific_marker_genes.R`: calculates dataset conserved marker genes across the subset of snRNA-seq datasets for cell types and subtypes

The scripts to generate the main and supplemental figures are included in the `figures/` directory.
- `plot_multispecies_integration.R`: plots the UMAP and marker gene plots for the main and supplemental figures
- `plot_multispecies_proportions.R`: calculates the and plots the cell type proportions across striatal subregions

![main integrated umaps of multi-omics, multi-species striatal single cell datasets](figures/main_fig6_multi-species_single_cell.jpg)
