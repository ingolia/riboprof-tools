#!/bin/bash

INBAM=NINM011_AGCTA_1M
GENOME=../data/sac_cer_yassour.bed

RUST_BACKTRACE=1 ~/Prog/riboprof-tools/target/debug/fp-framing \
            -l 17,41 \
            -b ${GENOME} \
            ~/Prog/riboprof-tools/test/${INBAM}.bam \
            -a ${INBAM}_rsannot.bam \
            -o ${INBAM}_rs
