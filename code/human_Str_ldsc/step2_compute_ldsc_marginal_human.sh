#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,short1,gpu
#SBATCH --time 2:00:00
#SBATCH --job-name=ldsc_mar
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/ldsc_mar_%A_%a.txt
#SBATCH --output=logs/ldsc_mar_%A_%a.txt
#SBATCH --array=1-122

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics/code/human_Str_ldsc 
DATADIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics/data/tidy_data/human_Str_ldsc
cd $CODEDIR; source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv )
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}//gwas_list_sumstats.tsv )
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv )

OUTDIR=${DATADIR}/enrichments; mkdir -p $OUTDIR

#################################################################################
# run LD score regression over the Stauffer striatum cell type binary annotations
CTS_FN1=${DATADIR}/Human_Striatum_snATAC.${POP}_hg38.ldcts
if [[ ! -f "$OUTDIR/Human_Striatum_snATAC.${GWAS_Label}.${POP}.cell_type_results.txt" ]]; then
ldsc.py --ref-ld-chr-cts $CTS_FN1 \
--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--h2-cts $GWAS --out $OUTDIR/Human_Striatum_snATAC.${GWAS_Label}.${POP}
fi
