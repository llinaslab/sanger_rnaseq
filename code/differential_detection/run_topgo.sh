#!/usr/bin/env bash

set -euo pipefail

go_terms=$1
desc=$2

# Lists comparing 3D7 and HB3
Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_3d7_hb3/off_3d7_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_3d7_hb3/off_3d7_not_hb3"

Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_3d7_hb3/off_hb3_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_3d7_hb3/off_hb3_not_3d7"

# Lists comparing 3D7 and IT
Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_3d7_it/off_3d7_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_3d7_it/off_3d7_not_it"

Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_3d7_it/off_it_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_3d7_it/off_it_not_3d7"

# Lists comparing HB3 and IT
Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_hb3_it/off_hb3_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_hb3_it/off_hb3_not_it"

Rscript code/topgo.R --gene_list "output/differential_detection/off_genes_hb3_it/off_it_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_genes_hb3_it/off_it_not_hb3"

# Genes only on in one strain
Rscript code/topgo.R --gene_list "data/goe/on_only_3d7_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/on_only_3d7"

Rscript code/topgo.R --gene_list "output/differential_detection/on_only_hb3_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/on_only_hb3"

Rscript code/topgo.R --gene_list "output/differential_detection/on_only_it_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/on_only_it"

# Genes only off in one strain
Rscript code/topgo.R --gene_list "output/differential_detection/off_only_3d7_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_only_3d7"

Rscript code/topgo.R --gene_list "output/differential_detection/off_only_hb3_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_only_hb3"

Rscript code/topgo.R --gene_list "output/differential_detection/off_only_it_list.txt" \
  --go_terms $go_terms \
	--anno $desc \
	--out_prefix "output/differential_detection/goe/off_only_it"

# Move all output plots to plots directory
#find data/goe -name *.pdf -exec cp -f {} plots \;
