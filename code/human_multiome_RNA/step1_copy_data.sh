PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/StriatumComparativeGenomics
DATADIR=${PROJDIR}/data/raw_data/STARsolo_out

PROJDIR2=/projects/pfenninggroup/singleCell/Logan_OUD-MDD-comorbid_ACC_NAC_multiome
DATADIR2=${PROJDIR2}/data/raw_data/STARsolo_out/all_runs

rsync -Pav $DATADIR2/*GE37* $DATADIR
rsync -Pav $DATADIR2/*GE25* $DATADIR
rsync -Pav $DATADIR2/*GE29* $DATADIR
rsync -Pav $DATADIR2/*GE41* $DATADIR
rsync -Pav $DATADIR2/*GE33* $DATADIR
