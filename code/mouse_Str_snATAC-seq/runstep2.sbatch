#!/bin/bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p pfen1
#SBATCH -o logs/step2_%j.txt
#SBATCH --mem=60GB

source activate ArchR_Env
Rscript step2_clustering_celltype.R