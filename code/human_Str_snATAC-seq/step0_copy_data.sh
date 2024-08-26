PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics
PROJDIR2=/projects/pfenninggroup/singleCell/BICCN_human_CATlas_snATAC-seq
PROJDIR3=/projects/pfenninggroup/singleCell/Corces2020_human_snATAC-seq

rsync -Paq ${PROJDIR2}/code/raw_code/human_Str_snATAC-seq \
${PROJDIR}/raw_code


rsync -Paq ${PROJDIR2}/data/raw_data/arrow/CaB* ${PROJDIR}/data/raw_data/arrow 
rsync -Paq ${PROJDIR2}/data/raw_data/arrow/NAC* ${PROJDIR}/data/raw_data/arrow 
rsync -Paq ${PROJDIR2}/data/raw_data/arrow/Pu.* ${PROJDIR}/data/raw_data/arrow 
rsync -Paq ${PROJDIR3}/data/arrow/*CAUD* ${PROJDIR}/data/raw_data/arrow 

