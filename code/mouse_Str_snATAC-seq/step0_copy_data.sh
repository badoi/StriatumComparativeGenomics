PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics
PROJDIR2=/projects/pfenninggroup/singleCell/BICCN_mouse_CATlas_snATAC-seq

rsync -Paq ${PROJDIR2}/data/raw_data/arrow/CP* ${PROJDIR}/data/raw_data/arrow_mm
rsync -Paq ${PROJDIR2}/data/raw_data/arrow/ACB* ${PROJDIR}/data/raw_data/arrow_mm

rsync -Paq ${PROJDIR2}/code/raw_code/mouse_Str_snATAC-seq_v2/* \
${PROJDIR}/code/mouse_Str_snATAC-seq

rsync -Paq ${PROJDIR2}/data/tidy_data/ArchRProjects/BICCN_mouse_Str_snATAC2 \
${PROJDIR}/data/tidy_data/ArchRProjects


