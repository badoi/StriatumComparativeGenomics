PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics
DATADIR=${PROJDIR}/data/tidy_data/mouse_merfish

PROJDIR2=/projects/pfenninggroup/singleCell/Zeng_transcriptome_Allen_10x_cells_wholebrainmouse_2024/data/tidy_data/merfish-objects-symbol

mkdir -p ${DATADIR}
rsync -Paq ${PROJDIR2}/*.spatial.seurat.rds ${DATADIR}
