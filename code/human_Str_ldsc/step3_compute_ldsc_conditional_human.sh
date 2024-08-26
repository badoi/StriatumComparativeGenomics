#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 1-0
#SBATCH --job-name=ldsc_con
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/ldsc_con_%A_%a.txt
#SBATCH --output=logs/ldsc_con_%A_%a.txt
#SBATCH --array=1-122

log2results() {
	awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
	# append line with leading white space to previous
	sed ':a;N;$!ba;s/\n //g'| \
	# parse Partitioned heritability rownames, starts w/ colon
	awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
	# transpose row to columns
	awk -v OFS='\t' '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
	    }
	}
	NF>p { p = NF }
	END {    
	    for(j=1; j<=p; j++) {
	        str=a[1,j]
	        for(i=2; i<=NR; i++){
	            str=str"\t"a[i,j];
	        }
	        print str
	    }
	}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
	gzip > ${OUT}.results.gz
}

# get the GWAS for this array job
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics
CODEDIR=${PROJDIR}/code/human_Str_ldsc 
DATADIR=${PROJDIR}/data/tidy_data/human_Str_ldsc
cd $CODEDIR; source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv )
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv )
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv )

OUTDIR=${DATADIR}/conditional; mkdir -p $OUTDIR

##########################################
## extract the annotations per cell type 
CTS_FN1=${DATADIR}/Human_Striatum_snATAC.${POP}_hg38.ldcts
CELLS=$( cut -f1 $CTS_FN1 | sed 's/Human_Striatum_snATAC.//g'| sed 's/.hg38//g' | sort | uniq )
BGPEAKS=$( cut -f2 $CTS_FN1 | cut -d ',' -f2 | sort | uniq )

## gather the foreground cell types for conditional cell type conditional
CELLTYPES=""
for CELL in $CELLS; do CELLTYPES=${CELLTYPES},$(awk -v VAR="${CELL}" '{if(match($1, VAR)) print $2}' $CTS_FN1 | cut -d ',' -f1 ) ; done
CELLTYPES=$(echo $CELLTYPES | sed 's/^,//g' | sed 's/ /,/g')

##############################################
## calculate the conditional gwas conditional
OUT=$OUTDIR/Human_Striatum_snATAC.conditional.${GWAS_Label}.${POP}
if [[ ! -f "${OUT}.results.gz" ]]; then 
ldsc.py \
--h2 $GWAS --print-coefficients \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${BGPEAKS},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
--out $OUT

log2results

fi


