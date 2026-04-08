#!/bin/bash

INBAM=NINM011_AGCTA_1M
GENOME=../data/sac_cer_yassour.bed
ASITES=../data/NINM011_asites.txt

RUST_BACKTRACE=1 ~/Prog/riboprof-tools/target/debug/fp-count \
            -b ${GENOME} \
            ~/Prog/riboprof-tools/test/${INBAM}_sorted.bam \
            -a ${ASITES} \
            -o ${INBAM}_rs_count.txt
