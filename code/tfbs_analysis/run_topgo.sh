#!/usr/bin/env bash

set -euo pipefail

go_terms=$1
desc=$2

# 3D7 activity
for file in $(ls output/tfbs_analysis/3d7_activity/*targets.txt); do
  Rscript code/topgo.R --gene_list "$file" \
    --go_terms $go_terms \
	  --anno $desc \
	  --out_prefix "${file%%.*}"
done

# HB3 activity
for file in $(ls output/tfbs_analysis/hb3_activity/*targets.txt); do
  Rscript code/topgo.R --gene_list "$file" \
    --go_terms $go_terms \
	  --anno $desc \
	  --out_prefix "${file%%.*}"
done

# IT activity
for file in $(ls output/tfbs_analysis/it_activity/*targets.txt); do
  Rscript code/topgo.R --gene_list "$file" \
    --go_terms $go_terms \
	  --anno $desc \
	  --out_prefix "${file%%.*}"
done
