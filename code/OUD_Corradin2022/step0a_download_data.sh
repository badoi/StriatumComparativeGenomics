#!/bin/bash

# download all roadmap data (use E070 for human dlpfc)
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
# unzip and clean
tar xzvf all.mnemonics.bedFiles.tgz
rm -r fall.mnemonics.bedFiles.tgz
# move files
mkdir -p data/roadmap
mv *.gz data/roadmap
# lift hg19 to hg38
$HOME/resources/./liftOver data/roadmap/E073_15_coreMarks_mnemonics.bed.gz $HOME/resources/hg19ToHg38.over.chain.gz data/roadmap/E073_15_coreMarks_mnemonics_hg38.bed data/roadmap/E073_15_coreMarks_mnemonics_hg38_unlifted.bed 
# write polycomb repressed regions
awk -F '\t' '$4 == "13_ReprPC" {print $1, $2, $3}' data/roadmap/E073_15_coreMarks_mnemonics_hg38.bed | sort -k 1,1 -k2,2n | sed 's/ /\t/g' > data/roadmap/E073_15_coreMarks_ReprPC_hg38_sorted.bed
