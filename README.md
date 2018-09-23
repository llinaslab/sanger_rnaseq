# Sanger RNA-seq

We wanted to identify UTRs in the falciparum genome using RNA-seq. A modified RNA-seq library prepartion protocol was applied to three strains of P. falciparum to predict this sites.

## ToDo

1. Figure 2 - Neighboring genes analysis
  - Calculate distances and correlations for all head-to-head, tail-to-tail, and head-to-tail genes for all neigboring protein-coding genes in the genome
  - Do this for translation translation start and stop sites and for newly predicted TSS start and stop sites.
  - For TSSa-RNAs as well?
  - Profile plots of example genes

2. Figure 3 - Predicted TSSs refine regulatory landscape
  - How many TSSs are dynamic, how many are secondary, how many are primary? How does that compare to the previous CAGE study done in Plasmodium? How well do our predictions line up with one another?
  - Do we see signs of sharp and/or broad promoters throughout the genome? What differentiates these kinds of promoters?
  - Are we able to use these UTR sequences to make better genome-wide predictions for active regulatory motifs?
  
3. Figure 5 - Comparison between the three strains
  - What transcripts can we detect in one strain, but in others? Two strains, but not in one? 
  - Are genes differentially expressed between strains?
  - Which genes are expressed at different time points throughout blood-stage development and why might that be?
  
## Workflow

### DAFT-seq

1. Align raw reads to Pf3D7 reference genome
2. Transcript discovery using bedtools and continuous blocks of coverage
3. Annotate discovered transcript (5' UTR, 3' UTR, ncRNA, etc.)
4. Transcript quantification using RSEM
5. Preprocessing for differential expression testing (low-count filter, bias removal, normalization)
6. Differential expression using EBSeq
7. Alternative splicing analysis using MISO + pre-defined splice sites based on alignment filters
8. Repeat for all strains, then compare results
9. Calculate phase for each transcript and compare between strains - which transcripts are out of phase from one another?
10. Calculate mean TPM and compare between strains - which transcripts are significantly more/less abundance in one strain over another?
11. Identify neighboring genes and the amount of intergenic space (or lack thereof) between them
12. Determine co-expression between neighboring genes in different orientations using Pearson correlation

### CAGE-seq

1. Map reads to Pf3D7 reference genome using SMALT (for soft-clipping feature)
2. Identify 5' most mapped base and count for each nucleotide, the number of 5' mapped bases
3. Normalize and analyze using CAGEr to identify TCs and PCs
4. Annotate TCs and PCs to identify which genes may have more than one TC and/or PC
5. Characterize promoter architecture such as width (shape) and nucleotide content in different areas of the genome

---

A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
