setwd("")
source("code/utils.R")

load_essentials()
sshhh("CAGEr")
sshhh("GenomicRanges")
sshhh("BSgenome.Pfalciparum.PlasmoDB.v24")

# collect command line args
args <- commandArgs(TRUE)

indir  <- args[1]
outdir <- args[2]

# input file paths
pathsToInputFiles <- list.files(indir, full.names = T, pattern = "*.ctss")

# create CAGEset object
cageset <- new("CAGEset",
               genomeName = "BSgenome.Pfalciparum.PlasmoDB.v24",
               inputFiles = pathsToInputFiles,
               inputFilesType = "ctss",
               sampleLabels = c("tp1", "tp2", "tp3", "tp4", "tp5", "tp6", "tp7"))

# read in files
getCTSS(cageset)

# change to outdir
setwd(outdir)

# plot raw correlation values
plotCorrelation(cageset, what = "CTSS", samples = "all", method = "pearson")

# plot reverse cumulative distributions per sample to estimate parameters to properly normalize
plotReverseCumulatives(cageset, values = "raw", onePlot = TRUE)

# normalize tag counts using the power law method
normalizeTagCount(cageset, method = "powerLaw",
                  fitInRange = c(10, 1000), alpha = 1.28, T = 1e5)

# plot reverse cumulative distributions per sample to see changes after normalization
plotReverseCumulatives(cageset,
                      values = "normalized", onePlot = TRUE)

# export the normalized counts to a bedgraph file
exportCTSStoBedGraph(cageset,
                    values = "normalized", oneFile = F)

# cluster TSSs into tag clusters (TCs)
clusterCTSS(object = cageset,
            threshold = 1,
            thresholdIsTpm = TRUE,
            nrPassThreshold = 1,
            method = "distclu",
            maxDist = 20,
            removeSingletons = TRUE,
            keepSingletonsAbove = 5,
            useMulticore = TRUE,
            nrCores = 8)

# calculate quantiles per tag cluster
cumulativeCTSSdistribution(cageset, clusters = "tagClusters")
quantilePositions(cageset, clusters = "tagClusters",
                  qLow = 0.1, qUp = 0.9)

# plot TC width
plotInterquantileWidth(cageset, clusters = "tagClusters",
                      tpmThreshold = 1, qLow = 0.1, qUp = 0.9)

# export to BED file
exportToBed(cageset, what = "tagClusters",
           qLow = 0.1, qUp = 0.9, oneFile = F)

# individual tag clusters
tc1 <- tagClusters(cageset, sample = "tp1",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc2 <- tagClusters(cageset, sample = "tp2",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc3 <- tagClusters(cageset, sample = "tp3",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc4 <- tagClusters(cageset, sample = "tp4",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc5 <- tagClusters(cageset, sample = "tp5",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc6 <- tagClusters(cageset, sample = "tp6",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
tc7 <- tagClusters(cageset, sample = "tp7",
                   returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)

# export tag clusters
export.gff3(GRanges(tc1), "tc1.gff")
export.gff3(GRanges(tc2), "tc2.gff")
export.gff3(GRanges(tc3), "tc3.gff")
export.gff3(GRanges(tc4), "tc4.gff")
export.gff3(GRanges(tc5), "tc5.gff")
export.gff3(GRanges(tc6), "tc5.gff")
export.gff3(GRanges(tc5), "tc5.gff")
export.gff3(GRanges(tc6), "tc6.gff")
export.gff3(GRanges(tc7), "tc7.gff")

# create consensus promoter clusters (PCs)
aggregateTagClusters(cageset,
                     tpmThreshold = 1,
                     qLow         = 0.1,
                     qUp          = 0.9,
                     maxDist      = 100,
                     excludeSignalBelowThreshold = F)

# calculate quantiles per PC
cumulativeCTSSdistribution(cageset, clusters = "consensusClusters")
quantilePositions(cageset, clusters = "consensusClusters",
                  qLow = 0.1, qUp = 0.9)

# plot PC width
#plotInterquantileWidth(cageset, clusters = "consensusClusters",
#                      tpmThreshold = 5, qLow = 0.1, qUp = 0.9)

# get expression profiles for consensus promoters
getExpressionProfiles(cageset, what = "consensusClusters",
                      tpmThreshold = 5, nrPassThreshold = 1, method = "som", xDim = 7, yDim = 1)
#plotExpressionProfiles(cageset, what = "consensusClusters")

# expression profile for indivdiual CTSSs
getExpressionProfiles(cageset, what = "CTSS",
                      tpmThreshold = 5, nrPassThreshold = 1, method = "som", xDim = 7, yDim = 1)
#plotExpressionProfiles(cageset, what = "CTSS")

# export to BED file with colors
exportToBed(cageset, what = "consensusClusters", colorByExpressionProfile = TRUE)
exportToBed(cageset, what = "CTSS", colorByExpressionProfile = TRUE)

# retrieve consensus clusters
cc <- consensusClusters(cageset)

cc$expression.cluster <- apply(cc, 1, function(x) {
  tmp <- str_replace_all(as.character(x[["consensus.cluster"]]), " ", "")
  if(tmp %in% names(cageset@consensusClustersExpressionClasses)) {
    cageset@consensusClustersExpressionClasses[[tmp]]
  } else {
    NA
  }
})

export.gff(cc, "promoter_clusters.gff")

#-------------------------------------------------------------------------
# Timecourse
#-------------------------------------------------------------------------

# Compare 1 to 2
scoreShift(cageset, groupX = "tp1", groupY = "tp2", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters12 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters12, "shifting_promoters12.tsv")

# Compare 2 to 3
scoreShift(cageset, groupX = "tp2", groupY = "tp3", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters23 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters23, "shifting_promoters23.tsv")

# Compare 3 to 4
scoreShift(cageset, groupX = "tp3", groupY = "tp4", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters34 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters34, "shifting_promoters34.tsv")

# Compare 4 to 5
scoreShift(cageset, groupX = "tp4", groupY = "tp5", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters45 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters45, "shifting_promoters45.tsv")

# Compare 5 to 6
scoreShift(cageset, groupX = "tp5", groupY = "tp6", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters56 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters56, "shifting_promoters56.tsv")

# Compare 6 to 7
scoreShift(cageset, groupX = "tp6", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters67 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters67, "shifting_promoters67.tsv")

#-------------------------------------------------------------------------
# Timepoint 1
#-------------------------------------------------------------------------

# Compare 1 to 3
scoreShift(cageset, groupX = "tp1", groupY = "tp3", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters13 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters13, "shifting_promoters13.tsv")

# Compare 1 to 4
scoreShift(cageset, groupX = "tp1", groupY = "tp4", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters14 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters14, "shifting_promoters14.tsv")

# Compare 1 to 5
scoreShift(cageset, groupX = "tp1", groupY = "tp5", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters15 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters15, "shifting_promoters15.tsv")

# Compare 1 to 6
scoreShift(cageset, groupX = "tp1", groupY = "tp6", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters16 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters16, "shifting_promoters16.tsv")

# Compare 1 to 7
scoreShift(cageset, groupX = "tp1", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters17 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters17, "shifting_promoters17.tsv")

#-------------------------------------------------------------------------
# Timepoint 2
#-------------------------------------------------------------------------

# Compare 2 to 4
scoreShift(cageset, groupX = "tp2", groupY = "tp4", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters24 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters24, "shifting_promoters24.tsv")

# Compare 2 to 5
scoreShift(cageset, groupX = "tp2", groupY = "tp5", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters25 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters25, "shifting_promoters25.tsv")

# Compare 2 to 6
scoreShift(cageset, groupX = "tp2", groupY = "tp6", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters26 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters26, "shifting_promoters26.tsv")

# Compare 2 to 7
scoreShift(cageset, groupX = "tp2", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters27 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters27, "shifting_promoters27.tsv")

#-------------------------------------------------------------------------
# Timepoint 3
#-------------------------------------------------------------------------

# Compare 3 to 5
scoreShift(cageset, groupX = "tp3", groupY = "tp5", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters35 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters35, "shifting_promoters35.tsv")

# Compare 3 to 6
scoreShift(cageset, groupX = "tp3", groupY = "tp6", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters36 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters36, "shifting_promoters36.tsv")

# Compare 3 to 7
scoreShift(cageset, groupX = "tp3", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters37 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters37, "shifting_promoters37.tsv")

#-------------------------------------------------------------------------
# Timepoint 4
#-------------------------------------------------------------------------

# Compare 4 to 6
scoreShift(cageset, groupX = "tp4", groupY = "tp6", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters46 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters46, "shifting_promoters46.tsv")

# Compare 4 to 7
scoreShift(cageset, groupX = "tp4", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters47 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters47, "shifting_promoters47.tsv")

#-------------------------------------------------------------------------
# Timepoint 5
#-------------------------------------------------------------------------

# Compare 5 to 7
scoreShift(cageset, groupX = "tp4", groupY = "tp7", testKS = TRUE, useTpmKS = TRUE)
shifting.promoters57 <- getShiftingPromoters(cageset, tpmThreshold = 5, fdrThreshold = 0.05)
write_tsv(shifting.promoters57, "shifting_promoters57.tsv")
