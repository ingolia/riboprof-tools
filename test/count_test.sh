#!/bin/bash

set -x
set -e

GENOME=../data/sac_cer_yassour.bed
ASITES=../data/NINM011_asites.txt

for BAM in `ls *_sort.bam`
do
    BASE=`basename ${BAM} .bam`

    RUST_BACKTRACE=1 ~/Prog/riboprof-tools/target/debug/fp-count \
		     -b ${GENOME} \
		     ~/Prog/riboprof-tools/test/${BASE}.bam \
		     -a ${ASITES} \
			 -d \
		     -o ${BASE}_rs_count.txt &
done

wait
