#!/usr/bin/env Rscript


# START -------------------------------------------------------------------

# clear environment
rm(list=ls())
options(echo=FALSE) # change to TRUE if you want arguments in output file

# quietly load packages
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}

# packages to load
auto_loads <- c("optparse",
                "readr",
                "tidyr",
                "dplyr",
                "cowplot",
                "stringr",
                "tools",
                "magrittr")

# load them quietly
invisible(lapply(auto_loads, sshhh))


# OPTIONS -----------------------------------------------------------------


options(warn=-1)

description <- "
This script will take in a normalized RPKM by transcript file and creates
plots based off this file to assess GC content bias. If you provide a list off
genes, with each gene separated by a new line character, then it will output
a PDF of each gene plotted instead of a boxplot
"

option_list <- list(
  make_option(c("-b", "--bias_file"), type="character",
              help=".rds or .tsv file", metavar="character"),
  make_option(c("-i", "--include"), type="character",
                help="List of genes to include", metavar="character"),
  make_option(c("-o", "--out_prefix"), type="character",
              help="Output file prefix", metavar="character"),
  make_option(c("-l", "--gene_list"), type="character",
              help="List of genes to plot", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,description=description)
opt <- parse_args(opt_parser)

if (is.null(opt$bias_file)) {
  print_help(opt_parser)
  stop("\nBias file is required\n", call.=FALSE)
}

if (is.null(opt$out_prefix)) {
  print_help(opt_parser)
  stop("\nOutput prefix is required is required\n", call.=FALSE)
}


# FUNCTIONS ---------------------------------------------------------------

include   <- read_tsv(opt$include,col_names=F)$X1
#gene_list <- ifelse(opt$gene_list,read_tsv(opt$gene_list,col_names=F)$X1,NA)

if (file_ext(opt$bias_file) == "tsv") {
  bias_file <- read_tsv(opt$bias_file,col_names=c("id","pos","tp","exp","gc"))
} else if (file_ext(opt$bias_file) == "rds") {
  bias_file <- readRDS(opt$bias_file)
} else {
  stop("Bias file must be either .tsv or .rds file. Aborting.")
}

# Plot gcbias scatter plots for individual transcripts
plot_gcbias_scatter <- function(df, gene) {

  df %>%
    filter(id == gene) %>%
    group_by(id, pos) %>%
    summarise(gc = unique(gc), exp = mean(exp)) %>%
    ggplot(aes(x=gc,y=exp)) +
    geom_point(alpha=0.75) +
    stat_smooth(color="red",method="lm") +
    ggtitle(gene) +
    ylab("Normalized Expression") +
    xlab("%GC")

}

# Plot gcbias scatter plots for multiple transcripts
plot_gcbias_scatter2 <- function(df) {

  df %>%
    group_by(id, pos) %>%
    summarise(gc=unique(gc),exp=mean(exp,na.rm=T)) %>%
    group_by(pos) %>%
    summarise(gc=mean(gc),exp=mean(exp,na.rm=T)) %>%
    ggplot(aes(x=gc,y=exp)) +
    geom_point(alpha=0.75) +
    stat_smooth(color="red",method="lm") +
    ggtitle(paste("Average of",length(unique(df$id)),"transcripts")) +
    ylab("Normalized Expression") +
    xlab("%GC")

}

# Plot gcbias across individual transcripts
plot_gcbias_transcript <- function(df, gene) {

  g1 <- df %>%
    filter(id == gene) %>%
    group_by(id, pos) %>%
    summarise(gc = unique(gc), exp = mean(exp)) %>%
    ggplot(aes(x=pos,y=exp)) +
    stat_smooth(color="blue",se=F,span=0.15) +
    ylab("Normalized Expression") +
    xlab("") +
    ggtitle(gene) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()) +
    geom_hline(yintercept = 1, linetype=2)

  g2 <- df %>%
    filter(id == gene) %>%
    group_by(id, pos) %>%
    summarise(gc = unique(gc), exp = mean(exp)) %>%
    ggplot(aes(x=pos,y=gc)) +
    stat_smooth(color="red",se=F,span=0.15) +
    ylab("%GC") +
    xlab("Position")

  plot_grid(g1,g2,nrow=2,align="v", rel_heights = c(2, 1))
}

# Plot gcbias across multiple transcripts
plot_gcbias_transcript2 <- function(df) {

  g1 <- df %>%
    group_by(id, pos) %>%
    summarise(exp=mean(exp,na.rm=T)) %>%
    group_by(pos) %>%
    summarise(exp=mean(exp,na.rm=T)) %>%
    ggplot(aes(x=pos,y=exp)) +
    stat_smooth(color="blue",se=F,span=0.25) +
    geom_point(color="blue") +
    ylab("Normalized Expression") +
    xlab("") +
    scale_x_continuous(limits=c(1,max(unique(df$pos)))) +
    ggtitle(paste("Average of",length(unique(df$id)),"transcripts")) +
    geom_hline(yintercept = 1, linetype=2)

  g2 <- df %>%
    group_by(id, pos) %>%
    summarise(gc=unique(gc)) %>%
    group_by(pos) %>%
    summarise(gc=mean(gc)) %>%
    ggplot(aes(x=pos,y=gc)) +
    stat_smooth(color="red",se=F,span=0.25) +
    geom_point(color="red") +
    ylab("%GC") +
    xlab("Position") +
    scale_x_continuous(limits=c(1,max(unique(df$pos))))

  plot_grid(g1,g2,nrow=2,align="v",rel_heights=c(1, 1))

}

# Plot gcbias across transcripts on top of each other
# Since ggplot2 doesn't let you create two plots using two different
# y-axes I needed to use the basic plot function in R
plot_gcbias_transcript3 <- function(df, gene) {

  if(gene %in% df$id) {

    what <- df %>%
      filter(id == gene) %>%
      spread(tp,exp)

    what$exp <- apply(what, 1, function(row) {
      mean(c(as.numeric(row[["tp1"]]),
             as.numeric(row[["tp2"]]),
             as.numeric(row[["tp3"]]),
             as.numeric(row[["tp4"]]),
             as.numeric(row[["tp5"]]),
             as.numeric(row[["tp6"]]),
             as.numeric(row[["tp7"]])))
    })

    if(sum(is.na(what$exp)) == 0) {

      par(mar = c(5,5,2,5))
      with(what, plot(pos, predict(loess(exp ~ pos,span=0.25)),
                      type="l",
                      col="blue",
                      ylab="Normalized Expression",
                      xlab="Transcript Position (5' - 3')",
                      lwd=2,
                      main=gene))
      abline(h=1,lty=2,col="grey2")

      par(new = T)
      with(what, plot(pos, predict(loess(gc ~ pos,span=0.25)),
                      type="l",col="red3",axes=F, xlab=NA, ylab=NA,lwd=1,lty=3))
      axis(side=4)
      mtext(side=4,line=3,"%GC")

    }

  } else {
    print(paste(gene, "NA"))
  }

}

# Bin the windows by GC content and plot boxplots of the normalized
# exprssion values within those bins
plot_gcbias_boxplot <- function(df, include) {

  # remove genes if only meant to include certain genes
  if(sum(!is.na(include))>0) {
    df <- df %>%
      filter(id %in% include)
  }

  df$bins <- as.character(cut(df$gc, breaks=c(0,0.1,0.2,0.3,0.4,0.5,1.0),
        labels=c("0-10%","11-20%","21-30%","31-40%","41-50%","51-100%")))

  df %>%
    select(-gc) %>%
    group_by(id, pos) %>%
    summarise(exp = mean(exp), bins = first(bins)) %>%
    ggplot(aes(x=bins,y=log2(exp))) +
    geom_boxplot() +
    xlab("%GC") +
    ylab("Log2 Normalized RPKM")

}


print_check <- function(what) {
  print(head(what))
  q(status=0)
}


# MAIN --------------------------------------------------------------------


# get transcripts we care about
tmp25 <- bias_file %>%
  group_by(id) %>%
  summarise(m = max(pos)) %>%
  filter(m < 25) %$%
  id

tmp50 <- bias_file %>%
  group_by(id) %>%
  summarise(m = max(pos)) %>%
  filter(m >= 25 & m <= 50) %$%
  id

tmp75 <- bias_file %>%
  group_by(id) %>%
  summarise(m = max(pos)) %>%
  filter(m >= 50 & m <= 75) %$%
  id

tmp100 <- bias_file %>%
  group_by(id) %>%
  summarise(m = max(pos)) %>%
  filter(m >= 75 & m <= 100) %$%
  id

tmpgt <- bias_file %>%
  group_by(id) %>%
  summarise(m = max(pos)) %>%
  filter(m > 100) %$%
  id

# some transcripts overlap...can't get position. Don't include these
ugh <- bias_file %>%
  group_by(id, pos) %>%
  summarise(l = length(unique(gc))) %>%
  filter(l > 1) %$%
  unique(id)

# plot gcbias as a scatter plot
g <- bias_file %>% filter(id %in% tmp25 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_scatter2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_scatter25.png"))
g <- bias_file %>% filter(id %in% tmp50 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_scatter2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_scatter50.png"))
g <- bias_file %>% filter(id %in% tmp75 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_scatter2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_scatter75.png"))
g <- bias_file %>% filter(id %in% tmp100 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_scatter2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_scatter100.png"))
g <- bias_file %>% filter(id %in% tmpgt & id %in% include & !(id %in% ugh)) %>% plot_gcbias_scatter2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_scattergt.png"))

# plot gcbias across transcript position
g <- bias_file %>% filter(id %in% tmp25 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_transcript2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_transcript25.png"))
g <- bias_file %>% filter(id %in% tmp50 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_transcript2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_transcript50.png"))
g <- bias_file %>% filter(id %in% tmp75 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_transcript2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_transcript75.png"))
g <- bias_file %>% filter(id %in% tmp100 & id %in% include & !(id %in% ugh)) %>% plot_gcbias_transcript2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_transcript100.png"))
g <- bias_file %>% filter(id %in% tmpgt & id %in% include & !(id %in% ugh)) %>% plot_gcbias_transcript2()
ggsave(plot = g, file = paste0(opt$out_prefix,"_transcriptgt.png"))

# Plot gcbias plot as aboxplot
g <- plot_gcbias_boxplot(bias_file,include)
ggsave(plot = g, file = paste0(opt$out_prefix,"_boxplot.png"))

q(status=0)
