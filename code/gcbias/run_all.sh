#!/usr/bin/env bash

set euo -pipefail

# ###############################################
# Create the 50bp window files for annotations in this directory using bedtools and awk
# ###############################################
bedtools makewindows -b ../data/gcbias/core_cds_exons_with_utrs.gff -w 50 | \
  bedtools intersect -a - -b ../data/gcbias/core_cds_exons_with_utrs.gff -wo | \
	awk '{print $1, "bedtools", "window", $2+1,$3,"50",$10,".",$12}' > ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.gff

bedtools makewindows -b ../data/gcbias/core_cds_exons.gff -w 50 | \
  bedtools intersect -a - -b ../data/gcbias/core_cds_exons.gff -wo | \
	awk '{print $1, "bedtools", "window", $2+1,$3,"50",$10,".",$12}' > ../data/gcbias/core_cds_exons_50bp_windows.gff

# ###############################################
# Get the nucleotide content for each window using bedtools
# ###############################################
bedtools nuc -fi ../data/gcbias/PlasmoDB-24_Pfalciparum3D7_Genome.fasta -bed ../data/gcbias/core_cds_exons_with_utrs.gff -s > ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.nuc
bedtools nuc -fi ../data/gcbias/PlasmoDB-24_Pfalciparum3D7_Genome.fasta -bed ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.gff -s > ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.nuc

# ###############################################
# Use bedtools to count reads over every interval per bam file
# ###############################################
bedtools multicov -bams ../bam/mapped_to_3d7/3d7.3d7_v3_chr.tp{1..7}.bam \
  -bed ../data/gcbias/core_cds_exons_50bp_windows.gff \
	-split -s -q 30 > ../data/gcbias/core_cds_exons_50bp_windows.counts

for strain in hb3 it; do
  bedtools multicov -bams ../bam/mapped_to_3d7/3d7.${strain}chr.tp{1..7}.bam \
	-bed ../data/gcbias/core_cds_exons_50bp_windows.gff \
	-split -s -q 30 > ../data/gcbias/${strain}_core_cds_exons_50bp_windows.counts;
done

bedtools multicov -bams ../bam/mapped_to_3d7/3d7.3d7_v3_chr.tp{1..7}.bam \
  -bed ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.gff \
	-split -s -q 30 > ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.counts

# ###############################################
# Add index to each window per transcript id
# ###############################################
awk '{if(n[$9] != 1) {n[$9] = 1; p = 1} else {p += 1} print $0,p}' ../data/gcbias/3d7_core_cds_exons_50bp_windows.counts > tmp && cat tmp > ../data/gcbias/3d7_core_cds_exons_50bp_windows.counts
awk '{if(n[$9] != 1) {n[$9] = 1; p = 1} else {p += 1} print $0,p}' ../data/gcbias/hb3_core_cds_exons_50bp_windows.counts > tmp && cat tmp > ../data/gcbias/hb3_core_cds_exons_50bp_windows.counts
awk '{if(n[$9] != 1) {n[$9] = 1; p = 1} else {p += 1} print $0,p}' ../data/gcbias/it_core_cds_exons_50bp_windows.counts > tmp && cat tmp > ../data/gcbias/it_core_cds_exons_50bp_windows.counts
awk '{if(n[$9] != 1) {n[$9] = 1; p = 1} else {p += 1} print $0,p}' ../data/gcbias/3d7_core_cds_exons_with_utrs_50bp_windows.counts > tmp && cat tmp > ../data/gcbias/3d7_core_cds_exons_with_utrs_50bp_windows.counts

# ###############################################
# Create processed data that's easy to plot and analyze
# ###############################################
for strain in 3d7 hb3 it; do Rscript ../scripts/gcbias/analyze_gcbias.R --counts ../data/gcbias/${strain}_core_cds_exons_50bp_windows.counts --nucs ../data/gcbias/core_cds_exons_50bp_windows.nuc; done
Rscript ../scripts/gcbias/analyze_gcbias.R --counts ../data/gcbias/3d7_core_cds_exons_with_utrs_50bp_windows.counts --nucs ../data/gcbias/core_cds_exons_with_utrs_50bp_windows.nuc

# ###############################################
# Create plots by runing the plot_gcbias.R script - plots are output to "plots/" directory
# ###############################################
Rscript ../code/gcbias/plot_gcbias.R --bias_file ../data/gcbias/3d7_core_cds_exons_with_utrs_50bp_windows.rds --include ../data/gcbias/on3d7.txt --out_prefix ../output/gcbias/3d7_with_utrs; done
for strain in hb3 it; do Rscript ../scripts/gcbias/plot_gcbias.R --bias_file ../data/gcbias/${strain}_core_cds_exons_50bp_windows.rds --include ../data/gcbias/on${strain}.txt --out_prefix output/gcbias/${strain}; done
