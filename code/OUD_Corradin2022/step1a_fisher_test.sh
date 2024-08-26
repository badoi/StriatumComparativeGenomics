#!/bin/bash

bedtools fisher -a data/NIHMS1778910-supplement-Table_S4_gain.bed -b data/roadmap/E073_15_coreMarks_ReprPC_hg38_sorted.bed -g ~/resources/hg38.chrom.sizes.sorted


bedtools fisher -a data/NIHMS1778910-supplement-Table_S4_lost.bed -b data/roadmap/E073_15_coreMarks_ReprPC_hg38_sorted.bed -g ~/resources/hg38.chrom.sizes.sorted
