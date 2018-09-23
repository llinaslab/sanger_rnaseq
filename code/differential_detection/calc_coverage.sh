#!/usr/bin/env bash

set -euo pipefail

if [[ $# < 2 ]]; then
echo ""
echo "Usage: bash $(basename $0) <num_threads> <annotation.gff>"
echo ""
echo "Run bedtools coverage on every alignment file for 3D7, HB3, and IT."
echo "This tool will invoke GNU paralell to calculate the coverage over every"
echo "interval within <annotations.gff>"
echo ""
exit 0
fi

threads=$1
anno=$2
bams=$3


parallel -j "${threads}" "bedtools coverage -a ${anno} -b ${bams}/3d7.3d7_v3_chr.tp{}.bam -s -split > ./data/cov/3d7_tp{}_cov.tsv" ::: {1..7}
parallel -j "${threads}" "bedtools coverage -a ${anno} -b ${bams}/{1}.3d7chr.tp{2}.bam -s -split > ./data/cov/{1}_tp{2}_cov.tsv" ::: hb3 it ::: {1..7}

exit 0
