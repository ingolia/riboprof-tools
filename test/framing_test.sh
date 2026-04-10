#!/bin/bash

set -x
set -e

GENOME=../data/sac_cer_yassour.bed

for BAM in `ls *_sort.bam`
do
    BASE=`basename ${BAM} .bam`
    RUST_BACKTRACE=1 ~/Prog/riboprof-tools/target/debug/fp-framing \
		     -l 17,41 \
		     -b ${GENOME} \
		     ~/Prog/riboprof-tools/test/${BASE}.bam \
		     -a ${BASE}_rsannot.bam \
		     -o ${BASE}_rs &

    RUST_BACKTRACE=1 ~/Prog/riboprof-tools/target/debug/fp-framing \
		     -l 27,29 \
		     -b ${GENOME} \
		     ~/Prog/riboprof-tools/test/${BASE}.bam \
		     -o ${BASE}_rs_narrow &

done

wait
