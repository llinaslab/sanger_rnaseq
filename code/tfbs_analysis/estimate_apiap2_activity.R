# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(NeatMap)
library(glmnet)
source("../code/tfbs_analysis/lasso.R")

main <- function(exp_file, mot_file, out_dir = "./", fdr = 0.01, cc = 0.5) {

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  # read in expression and fitted data
  e <- as.data.frame(read.csv(exp_file, sep="\t", row.names=1, header=T))
  m <- as.matrix(read.csv(mot_file, sep="\t", row.names=1, header=T))
  # find activity profiles
  coefs <- PBM.glmnet.getCoefs(e, m, file = paste0(out_dir, "/AP2-coefs.txt"))
  # sort activity profiles
  PBM.drawCoefHeatMap(coefs, file = paste0(out_dir, "/AP2-coefs.txt"))
  # ...what?
  PBM.glmnet.getFit(e, m, file=paste0(out_dir, "/fit-cores-enet-weighted.eps"))
  # fdr=0.01 is stringent; in many cases,
  # eg if you have few expression timepoint,
  # you can relax it down to fdr=0.2 (or higher but be careful of potential FP)
  PBM.writeTargetFiles(coefs, e, m, fdr=fdr, corr=cc, outdir=out_dir)
  PBM.plotAP2CoefVSExp(coefs, e, m, fdr=fdr, corr=cc, file=paste0(out_dir,"/coef-vs-exp.pdf"))

}

dir <- "../output/tfbs_analysis/"
strains <- c("3d7","hb3", "it")
false_discovery_rate <- 0.2
correlation_cutoff <- 0.5

for (strain in strains) {
  main(paste0(dir, strain, "_tpms.exp"),
       paste0(dir, strain, "_scores.mot"),
       out_dir = paste0(dir, strain, "_activity"),
       fdr = false_discovery_rate,
       cc = correlation_cutoff)
}

notheme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),legend.position="none",
                 panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),plot.background=element_blank())

color <- scale_fill_gradient2(low="blue", high="yellow", mid="black", midpoint=0)

for (strain in c("3d7","hb3","it")) {

filePath <- paste0("../output/tfbs_analysis/",strain,"_activity/AP2-coefs.txt")
outDir <- paste0("../output/tfbs_analysis/",strain,"_activity/")

m <- read.table(filePath, sep="\t", stringsAsFactors = F, header = T)
row.names(m) <- m$X
m$X <- NULL
colnames(m) <- c("T0", "T8", "T16", "T24", "T32", "T40", "T48")

mds <- nMDS(m)$x
o <- as.data.frame(apply(mds[,1:2], 1, function(x) {atan2(x[1], x[2])}))
colnames(o) <- "exp"
o$motif <- row.names(o)
row.names(o) <- NULL
o <- o[with(o, order(o$exp)),]$motif

m$motif <- row.names(m)
row.names(m) <- NULL
m %<>% gather(tp, exp, -motif)

m <- transform(m, motif = factor(motif, levels = rev(o)))

g <- ggplot(m, aes(x = tp, y = motif)) +
  geom_tile(aes(fill=exp)) +
  scale_fill_gradient2(low="blue", mid="black", high="yellow", space="rgb", midpoint = 0, limits=c(-2,2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Time") +
  ylab("Motif") +
  theme_classic() +
  ggtitle(paste0("Promoter Cov ApiAP2 Activity"))

ggsave(plot = g, file = paste0(outDir, "AP2-coefs_hm.png"))

}
